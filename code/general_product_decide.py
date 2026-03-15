"""Compatibility shim for the relocated general_product_decide module.

Execution from this legacy path is forwarded to src/active/general_product_decide.py.
Imports from this legacy path are also forwarded so existing tooling keeps working.
"""

from __future__ import annotations

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import runpy

_TARGET = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "active"
    / "general_product_decide.py"
)


def _load_target_module():
    spec = spec_from_file_location("src.active.general_product_decide", _TARGET)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load target module: {_TARGET}")
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


if __name__ == "__main__":
    runpy.run_path(str(_TARGET), run_name="__main__")
else:
    _module = _load_target_module()
    for _name, _value in _module.__dict__.items():
        if _name.startswith("__") and _name not in {"__doc__", "__all__"}:
            continue
        globals()[_name] = _value
