"""Convert ImageJ ROI .zip files into the _ROIs.npz format consumed by the
Soma Detector (v2.1 and v3.0).

Usage:  python convert_roi_zip_to_npz.py [folder]
The folder holds the *_ROI.zip files next to their .lsm / .czi stacks.  Without
an argument a folder-picker dialog opens (or set INPUT_DIR below).

For every <basename>_ROI.zip in the input directory:
  1. Locate the matching <basename>.lsm (or .czi) and read its stack shape
     and pixel size.
  2. Read each ROI; parse z-slice from the ROI name ("ZZZZ-YYYY-XXXX");
     rasterize the polygon/freehand/oval perimeter into a 2D boolean mask.
  3. Compute centroid (cy, cx) in pixels, area in µm².
  4. Pack everything into <basename>_ROIs.npz with the keys the soma
     detector expects (ch_idx, z, z_min, z_max, cy, cx, area_um2,
     cell_mean, peak_intensity, mask_<i>, rejected, source, stack_shape,
     px_um, czi_fname).

By default ROIs are assigned to channel 0 (GFP-Cre / green channel).
Intensities (cell_mean, peak_intensity) are filled with zeros — they get
re-computed when the soma detector exports measurements.
"""
import os
import re
import sys
import glob
import xml.etree.ElementTree as ET
import numpy as np
import roifile
import tifffile
import czifile
from skimage.draw import polygon as sk_polygon

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

INPUT_DIR = None             # optional fixed folder; None = command line / dialog
SKIP_IF_EXISTS = True       # don't overwrite an existing _ROIs.npz
DEFAULT_CHANNEL = 0          # ROIs land on green/GFP-Cre channel by default
DEFAULT_PX_UM = 0.4457       # fallback if pixel size cannot be read

# ── LSM metadata helpers ──────────────────────────────────────────────────────

def read_czi_meta(czi_path):
    """Return (stack_shape (c, z, y, x), px_um) for a Zeiss .czi file."""
    try:
        with czifile.CziFile(czi_path) as czi:
            arr = czi.asarray()
            meta_xml = czi.metadata()
        # Pixel size from CZI XML (Scaling/Items/Distance with Id="X")
        px_um = DEFAULT_PX_UM
        try:
            root = ET.fromstring(meta_xml)
            for d in root.iter('Distance'):
                if d.get('Id') == 'X':
                    val = d.find('Value')
                    if val is not None and val.text:
                        # CZI stores in metres
                        px_um = float(val.text) * 1e6
                        break
        except Exception:
            pass

        # CZI .asarray() returns shape like (1, C, Z, Y, X, 1) or similar.
        # Reduce singleton dims, then map to (C, Z, Y, X).
        a = np.squeeze(arr)
        if a.ndim == 4:
            stack_shape = tuple(int(x) for x in a.shape)  # assume (C, Z, Y, X)
        elif a.ndim == 3:
            # Could be (Z, Y, X) single channel or (C, Y, X) single slice — assume single C
            stack_shape = (1, int(a.shape[0]), int(a.shape[1]), int(a.shape[2]))
        elif a.ndim == 2:
            stack_shape = (1, 1, int(a.shape[0]), int(a.shape[1]))
        else:
            stack_shape = tuple(int(x) for x in a.shape)
        return stack_shape, px_um
    except Exception as e:
        sys.stderr.write(f"  CZI read error ({os.path.basename(czi_path)}): {e}\n")
        return None, DEFAULT_PX_UM


def read_lsm_meta(lsm_path):
    """Return (stack_shape (c, z, y, x), px_um). Falls back to (None, default)
    if the file cannot be opened."""
    try:
        with tifffile.TiffFile(lsm_path) as tf:
            # LSM stores a single series; .asarray() shapes it (Z, C, Y, X)
            # or (C, Z, Y, X) depending on writer. Inspect the LSM metadata.
            arr = tf.asarray()
            lsm_meta = tf.lsm_metadata or {}
            # Pixel size — values are in metres in the LSM header
            x_um = float(lsm_meta.get("VoxelSizeX", 0)) * 1e6
            y_um = float(lsm_meta.get("VoxelSizeY", 0)) * 1e6
            px_um = x_um if x_um > 0 else (y_um if y_um > 0 else DEFAULT_PX_UM)

            # Reorder dimensions to (C, Z, Y, X) — same convention the soma
            # detector uses internally.
            shape = arr.shape
            if len(shape) == 4:
                # Common variants: (Z, C, Y, X) or (C, Z, Y, X)
                # The LSM metadata DimensionChannels / DimensionZ tell us which is which.
                nC = int(lsm_meta.get("DimensionChannels", 1))
                nZ = int(lsm_meta.get("DimensionZ", 1))
                if shape[0] == nC and shape[1] == nZ:
                    stack_shape = shape
                elif shape[0] == nZ and shape[1] == nC:
                    stack_shape = (shape[1], shape[0], shape[2], shape[3])
                else:
                    # Best guess: assume (Z, C, Y, X) since that's what
                    # tifffile typically returns for LSM
                    stack_shape = (shape[1], shape[0], shape[2], shape[3])
            elif len(shape) == 3:
                # (C, Y, X) — single Z slice
                stack_shape = (shape[0], 1, shape[1], shape[2])
            elif len(shape) == 2:
                # (Y, X) — single channel, single Z
                stack_shape = (1, 1, shape[0], shape[1])
            else:
                stack_shape = shape

            return tuple(int(s) for s in stack_shape), px_um
    except Exception as e:
        sys.stderr.write(f"  LSM read error ({os.path.basename(lsm_path)}): {e}\n")
        return None, DEFAULT_PX_UM


# ── ROI → mask conversion ─────────────────────────────────────────────────────

def parse_z_from_name(name, default_z=0):
    """ImageJ multi-slice ROI names are typically "ZZZZ-YYYY-XXXX".
    Return the leading integer; ImageJ uses 1-based slice numbering, so
    convert to 0-based by subtracting 1."""
    parts = str(name).split("-")
    try:
        return max(0, int(parts[0]) - 1)
    except (ValueError, IndexError):
        return default_z


def roi_to_mask(roi, image_height, image_width):
    """Rasterize an ImageJ ROI to a (H, W) bool mask plus return its
    centroid and area-in-pixels."""
    coords = roi.coordinates()
    if coords is None or len(coords) < 3:
        # Fall back to bounding rectangle
        mask = np.zeros((image_height, image_width), dtype=bool)
        t, l, b, r = int(roi.top), int(roi.left), int(roi.bottom), int(roi.right)
        b = min(b, image_height - 1); r = min(r, image_width - 1)
        if b > t and r > l:
            mask[t:b + 1, l:r + 1] = True
        ys, xs = np.where(mask)
        if len(ys) == 0:
            return mask, float(t), float(l), 0
        return mask, float(ys.mean()), float(xs.mean()), int(mask.sum())

    # coords is (N, 2) of (x, y) — note ImageJ stores X first, Y second
    xs = np.asarray(coords[:, 0], dtype=np.float64)
    ys = np.asarray(coords[:, 1], dtype=np.float64)

    rr, cc = sk_polygon(ys, xs, shape=(image_height, image_width))
    mask = np.zeros((image_height, image_width), dtype=bool)
    mask[rr, cc] = True

    if mask.sum() == 0:
        # Polygon collapsed (single point / line) — fall back to bbox
        t, l, b, r = int(roi.top), int(roi.left), int(roi.bottom), int(roi.right)
        b = min(max(b, t + 1), image_height - 1)
        r = min(max(r, l + 1), image_width - 1)
        mask[t:b + 1, l:r + 1] = True

    ys_in, xs_in = np.where(mask)
    return mask, float(ys_in.mean()), float(xs_in.mean()), int(mask.sum())


# ── Conversion driver ────────────────────────────────────────────────────────

def find_image_for(roi_zip):
    """Find the matching image file by stripping the '_ROI' / '_roi' suffix
    (case-insensitive) and trying common extensions."""
    base, _ = os.path.splitext(roi_zip)
    base = re.sub(r"_ROI$", "", base, flags=re.IGNORECASE)
    for ext in (".lsm", ".czi", ".tif", ".tiff"):
        candidate = base + ext
        if os.path.isfile(candidate):
            return candidate
    return None


def read_image_meta(image_path):
    """Dispatch to the appropriate reader based on extension."""
    ext = os.path.splitext(image_path)[1].lower()
    if ext == ".czi":
        return read_czi_meta(image_path)
    return read_lsm_meta(image_path)


def convert_one(roi_zip, output_dir=None):
    """Convert one ROI.zip → _ROIs.npz. Returns the output path or None."""
    if output_dir is None:
        output_dir = os.path.dirname(roi_zip)

    base, _ = os.path.splitext(os.path.basename(roi_zip))
    base = re.sub(r"_ROI$", "", base, flags=re.IGNORECASE)
    out_path = os.path.join(output_dir, base + "_ROIs.npz")
    if SKIP_IF_EXISTS and os.path.isfile(out_path):
        return f"SKIP (exists): {os.path.basename(out_path)}"

    # Pair with the image file
    image_path = find_image_for(roi_zip)
    if image_path is None:
        return f"ERR  no matching image for {os.path.basename(roi_zip)}"

    stack_shape, px_um = read_image_meta(image_path)
    if stack_shape is None:
        return f"ERR  could not read LSM {os.path.basename(image_path)}"

    _n_ch, n_z, h, w = stack_shape

    rois = roifile.ImagejRoi.fromfile(roi_zip)
    if not rois:
        return f"EMPTY {os.path.basename(roi_zip)} contains 0 ROIs"

    # Build the per-ROI arrays
    ch_idx_arr  = []
    z_arr       = []
    cy_arr      = []
    cx_arr      = []
    area_arr    = []
    cell_mean   = []
    peak_int    = []
    rejected    = []
    source      = []
    masks       = []

    pixel_area_um2 = px_um * px_um

    for roi in rois:
        z = parse_z_from_name(roi.name, default_z=0)
        z = min(max(z, 0), n_z - 1)  # clip into stack range
        mask, cy, cx, area_px = roi_to_mask(roi, h, w)
        if area_px == 0:
            continue
        area_um2 = area_px * pixel_area_um2
        ch_idx_arr.append(DEFAULT_CHANNEL)
        z_arr.append(z)
        cy_arr.append(cy)
        cx_arr.append(cx)
        area_arr.append(area_um2)
        cell_mean.append(0.0)
        peak_int.append(0.0)
        rejected.append(False)
        source.append("manual")  # came from ImageJ user
        masks.append(mask)

    if not masks:
        return f"EMPTY all ROIs in {os.path.basename(roi_zip)} collapsed to 0 area"

    # Pack into the soma-detector npz layout
    n_rois = len(masks)
    save_kwargs = {
        "ch_idx":         np.asarray(ch_idx_arr, dtype=np.int32),
        "z":              np.asarray(z_arr, dtype=np.int32),
        "z_min":          np.asarray(z_arr, dtype=np.int32),
        "z_max":          np.asarray(z_arr, dtype=np.int32),
        "cy":             np.asarray(cy_arr, dtype=np.float64),
        "cx":             np.asarray(cx_arr, dtype=np.float64),
        "area_um2":       np.asarray(area_arr, dtype=np.float64),
        "cell_mean":      np.asarray(cell_mean, dtype=np.float64),
        "peak_intensity": np.asarray(peak_int, dtype=np.float64),
        "rejected":       np.asarray(rejected, dtype=bool),
        "source":         np.asarray(source, dtype="<U10"),
        "stack_shape":    np.asarray(stack_shape, dtype=np.int32),
        "px_um":          np.asarray(px_um, dtype=np.float64),
        "czi_fname":      np.asarray(os.path.basename(image_path), dtype="<U255"),
    }
    for i, m in enumerate(masks):
        save_kwargs[f"mask_{i}"] = m  # bool 2D

    np.savez_compressed(out_path, **save_kwargs)
    return f"OK   {n_rois:3d} ROIs -> {os.path.basename(out_path)}  (stack {stack_shape}, px {px_um:.4f}µm)"


def resolve_input_dir(argv):
    """Folder holding the *_ROI.zip files: the first command-line argument,
    else INPUT_DIR if set, else a folder-picker dialog (the current working
    directory when no display is available)."""
    if len(argv) > 1:
        return argv[1]
    if INPUT_DIR:
        return INPUT_DIR
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
    except Exception:
        return os.getcwd()
    chosen = filedialog.askdirectory(title="Folder with the ImageJ *_ROI.zip files")
    root.destroy()
    if not chosen:
        sys.exit("No folder chosen.")
    return chosen


def main():
    input_dir = resolve_input_dir(sys.argv)
    if not os.path.isdir(input_dir):
        sys.exit(f"Not a folder: {input_dir}")
    roi_zips = sorted(glob.glob(os.path.join(input_dir, "*_ROI.zip")))
    print(f"Found {len(roi_zips)} ROI.zip files in {input_dir}\n")
    for zfile in roi_zips:
        result = convert_one(zfile)
        print(result)
    print("\nDone.")


if __name__ == "__main__":
    main()
