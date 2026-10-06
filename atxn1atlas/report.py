"""Pretty-printing helpers."""


def table(rows, cols, title=None, pct=()):
    if title:
        print(f"\n{title}")
        print("-" * len(title))
    widths = {c: max(len(str(c)), *(len(_fmt(r.get(c), c in pct)) for r in rows)) for c in cols}
    print("  ".join(str(c).ljust(widths[c]) for c in cols))
    for r in rows:
        print("  ".join(_fmt(r.get(c), c in pct).ljust(widths[c]) for c in cols))


def _fmt(v, as_pct):
    if v is None:
        return ""
    if as_pct and isinstance(v, float):
        return f"{v * 100:.2f}%"
    if isinstance(v, float):
        return f"{v:.4f}"
    if isinstance(v, int):
        return f"{v:,}"
    return str(v)


def header(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)
