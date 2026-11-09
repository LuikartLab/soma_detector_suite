"""
cell_jump_ipc.py
----------------
Lightweight IPC for jumping between PTEN-FMR analysis GUIs.

The Soma Detector GUI (soma_area_and_intensity_gui_2.1.py) starts a
JumpServer at launch.  The Group Analysis GUI sends a short JSON command
via send_jump() when the user clicks "Open in soma detector".

Protocol
========
* TCP, localhost only (`127.0.0.1` by default).
* Each connection carries one newline-terminated JSON object.
* No reply is sent — fire-and-forget.

Message schema
==============
    {
        "action":      "open_cell",
        "file":        "<basename>.czi",   # required
        "search_dirs": [<absolute paths>], # optional list of dirs to search
        "slice":       <int>,              # optional, defaults to centred
        "soma_num":    <int>,              # informational
        "centroid_x_um": <float>,          # optional — pixel-matching anchor
        "centroid_y_um": <float>,
        "channel_hint": <int>,             # optional — which channel to display
        "condition":     "<str>"           # informational only
    }

Use 127.0.0.1 (loopback) so Windows Firewall does not flag the listener
and so no other host can reach it.  If port 5099 is in use the server
fails quietly and the group GUI's send_jump() returns False — its
fallback is then to copy the identifier to the clipboard.
"""
import socket
import threading
import json
import queue


DEFAULT_PORT = 5099
DEFAULT_HOST = '127.0.0.1'


# ── Client side (group analysis GUI) ──────────────────────────────────────────

def send_jump(message, port=DEFAULT_PORT, host=DEFAULT_HOST, timeout=1.0):
    """Send a single jump command to the soma GUI's listener.

    Returns True if the message was delivered (server accepted the
    connection and read the bytes), False otherwise.  Callers should treat
    False as "soma GUI is probably not running" and fall back to copying
    the identifier to the clipboard.
    """
    try:
        with socket.create_connection((host, port), timeout=timeout) as s:
            payload = (json.dumps(message) + '\n').encode('utf-8')
            s.sendall(payload)
        return True
    except OSError:
        return False


# ── Server side (soma detector GUI) ──────────────────────────────────────────

class JumpServer:
    """Background-thread TCP listener that hands messages to a Tk root via
    a thread-safe queue.

    The listener thread only does socket I/O and JSON parsing.  All Tk
    interaction happens on the main thread by way of root.after(), keeping
    the GUI thread-safe.

    Usage::

        srv = JumpServer(root, on_message=app.handle_jump)
        ok, err = srv.start()
        if not ok: print('IPC unavailable:', err)
        # … run mainloop …
        srv.stop()
    """

    def __init__(self, tk_root, on_message,
                 port=DEFAULT_PORT, host=DEFAULT_HOST):
        self.tk_root    = tk_root
        self.on_message = on_message
        self.port       = port
        self.host       = host
        self._q         = queue.Queue()
        self._stop      = threading.Event()
        self._sock      = None
        self._thread    = None
        self._drain_id  = None

    # ── public ──

    def start(self):
        """Bind, listen, and start the accept loop in a daemon thread.

        Returns (success: bool, error_string: str).
        """
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            s.bind((self.host, self.port))
        except OSError as e:
            s.close()
            return False, str(e)
        s.listen(2)
        s.settimeout(0.5)   # so the accept loop can check _stop periodically
        self._sock = s
        self._thread = threading.Thread(target=self._accept_loop, daemon=True)
        self._thread.start()
        self._schedule_drain()
        return True, ''

    def stop(self):
        """Stop the listener and close the socket.  Idempotent."""
        self._stop.set()
        try:
            if self._sock is not None:
                self._sock.close()
        except OSError:
            pass
        if self._drain_id is not None:
            try:
                self.tk_root.after_cancel(self._drain_id)
            except Exception:
                pass
            self._drain_id = None

    # ── internals ──

    def _accept_loop(self):
        """Run in a daemon thread.  Accepts one connection at a time and
        reads up to a newline.  Parsed messages are pushed to the queue;
        bad messages are silently dropped."""
        while not self._stop.is_set():
            try:
                conn, _addr = self._sock.accept()
            except (socket.timeout, OSError):
                continue
            try:
                conn.settimeout(2.0)
                data = b''
                while b'\n' not in data and len(data) < 65536:
                    chunk = conn.recv(4096)
                    if not chunk:
                        break
                    data += chunk
                if data.strip():
                    msg = json.loads(data.decode('utf-8').strip())
                    self._q.put(msg)
            except (OSError, ValueError, UnicodeDecodeError):
                pass
            finally:
                try:
                    conn.close()
                except OSError:
                    pass

    def _schedule_drain(self):
        """Drain the queue on the Tk main thread and re-schedule."""
        try:
            while True:
                msg = self._q.get_nowait()
                try:
                    self.on_message(msg)
                except Exception:
                    import traceback
                    traceback.print_exc()
        except queue.Empty:
            pass
        if not self._stop.is_set():
            self._drain_id = self.tk_root.after(150, self._schedule_drain)
