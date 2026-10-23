#!/usr/bin/env python3
"""HTML report of packbench results: packbench_report.py <results-dir> <out.html>

Reads <dir>/<d>_build.tsv (packbench build / hash output), <dir>/quick_<d>.tsv (query rows, each
prefixed with the dataset name) and, if present, <dir>/strategies_d2.tsv."""
import html
import os
import sys

DATASETS = {
    "d1": ("D1 · 450 mixed bacterial genomes", "2.05 Gbp in 66,010 contigs: many species, few copies of a k-mer"),
    "d2": ("D2 · 300 Bacteroides uniformis genomes", "1.31 Gbp in 45,328 contigs, one species: a k-mer occurs in many genomes (~150 seeds per lookup)"),
    "d3": ("D3 · human GRCh37 chr1–8", "1.54 Gbp in 8 chromosomes: repeats"),
}
# (kind, mode, thp, label, css color or None)
VARIANTS = [
    ("old", "heap", "0", "Flexmap blob (flexalign today)", "var(--old)"),
    ("old", "heap", "1", "Flexmap blob + huge pages", None),
    ("old", "mmap", "1", "Flexmap blob, mmap", None),
    ("compact", "heap", "0", "Packed, compact keys, no huge pages", None),
    ("compact", "heap", "1", "Packed, compact keys", "var(--new)"),
    ("compact", "mmap", "1", "Packed, compact keys, mmap", None),
    ("hashed", "heap", "1", "Packed, hashed keys", "var(--hash)"),
]
BASE = ("old", "heap", "0")
PACKED = [("compact", "heap", "1"), ("hashed", "heap", "1")]
# (set, work, label, column head)
WORKS = [
    ("perread", "seeds", "Short reads (150 bp), one call per read", "150 bp<br>per read"),
    ("long", "seeds", "Long reads (10 kb), one call per read", "10 kb<br>per read"),
    ("reads", "seeds", "Short reads as one stream", "150 bp<br>stream"),
    ("reads", "values", "All values of each key", "all values"),
    ("random", "values", "Random keys (mostly misses)", "random<br>keys"),
]


def rows(path):
    if not os.path.exists(path):
        return []
    return [l.rstrip("\n").split("\t") for l in open(path) if l.count("\t") > 1]


def main():
    rdir, out = sys.argv[1], sys.argv[2]
    builds, loads, times = {}, {}, {}
    for d in DATASETS:
        b = {}
        for f in rows(os.path.join(rdir, f"{d}_build.tsv")):
            b.setdefault(f[0], {}).update(dict(zip(f[2::2], f[3::2])) if f[0] == "build" else dict(zip(f[1::2], f[2::2])))
        if b:
            builds[d] = b
        for f in rows(os.path.join(rdir, f"quick_{d}.tsv")):
            key = (f[0], f[2], f[3], f[4].split("=")[1])
            if f[1] == "load":
                loads[key] = {"load_s": float(f[8]), "rss_gb": float(f[10])}
            elif f[1] == "time":
                times[key + (f[5], f[6], int(f[8]), int(f[10]))] = {"q": int(f[12]), "mqs": float(f[18]), "values": int(f[22]), "chk": f[24]}
    present = [d for d in DATASETS if any(k[0] == d for k in loads)]
    threads = sorted({k[6] for k in times})

    def t(d, v, w, th, pf=1):
        return times.get((d,) + tuple(v) + tuple(w[:2]) + (th, pf))

    def mismatches(d):
        bad = []
        for w in WORKS:
            for th in threads:
                chks = {x["chk"] for v in VARIANTS for pf in (0, 1) if (x := t(d, v[:3], w, th, pf))}
                if len(chks) > 1:
                    bad.append(f"{w[2]} at {th} threads")
        return bad

    h = []
    a = h.append
    a('<title>Packed Flexmap Benchmark</title>\n<meta name="claude:viewport-fit-sides" content="cover">')
    a('<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans+Condensed:wght@600&family=IBM+Plex+Sans:wght@400;600&display=swap">')
    a(f"<style>{CSS}</style><div class=\"wrap\">")
    a('<h1>Packed Flexmap Benchmark</h1><p class="lede">flexmap branch <span class="mono">packed-index</span> against <span class="mono">origin/zeroable_load</span> (cebfd56, the rev flexalign pins). '
      'The packed index holds the same k-mer positions and flank cells in less memory and answers every lookup with the same seeds. '
      'Quick runs: 2 M queries per dataset and scenario, one timed round after a warm-up.</p>')

    a('<h2>At a glance</h2><p class="muted">Short reads, one lookup call per 150 bp read, M seed lookups per second, both sides with their batch prefetch. Memory is the whole index in RAM.</p>')
    a('<div class="scroll"><table><tr><th>Dataset</th><th class="num">Old index</th><th class="num">Packed compact</th><th class="num">Packed hashed</th>')
    for th in threads:
        a(f'<th class="num">Old {th}t</th><th class="num">Compact {th}t</th><th class="num">Hashed {th}t</th>')
    a('<th>Same answers</th></tr>')
    for d in present:
        b = builds.get(d, {})

        def mem(k):
            return f'{float(b[k]["total_gb"]):.2f} GB' if k in b and "total_gb" in b[k] else "–"

        def pct(k):
            if k in b and "total_gb" in b[k] and "old" in b:
                return f' <span class="up">−{(1 - float(b[k]["total_gb"]) / float(b["old"]["total_gb"])) * 100:.0f}%</span>'
            return ""
        a(f'<tr><td>{html.escape(DATASETS[d][0])}</td><td class="num">{mem("old")}</td><td class="num">{mem("compact")}{pct("compact")}</td><td class="num">{mem("hashed")}{pct("hashed")}</td>')
        for th in threads:
            x0 = t(d, BASE, WORKS[0], th)
            a(f'<td class="num">{x0["mqs"]:.2f}</td>' if x0 else '<td class="num">–</td>')
            for v in PACKED:
                x = t(d, v, WORKS[0], th)
                a(f'<td class="num">{x["mqs"]:.2f} <span class="up">{x["mqs"] / x0["mqs"]:.2f}×</span></td>' if x and x0 else '<td class="num">–</td>')
        a(f'<td>{"<span class=ok>yes</span>" if not mismatches(d) else "<span class=fail>no</span>"}</td></tr>')
    a("</table></div>")

    for d in present:
        name, desc = DATASETS[d]
        a(f'<h2>{html.escape(name)}</h2><p class="muted">{desc}</p>')
        b = builds.get(d, {})
        tiles = []
        if "old" in b:
            parts = [f'old: keys {float(b["old"]["keys_gb"]):.2f} + values {float(b["old"]["values_gb"]):.2f}']
            v = f'{float(b["old"]["total_gb"]):.2f}'
            for k in ("compact", "hashed"):
                if k in b and "total_gb" in b[k]:
                    v += f' → {float(b[k]["total_gb"]):.2f}'
                    parts.append(f'{k}: keys {float(b[k]["keys_gb"]):.2f} + values {float(b[k]["values_gb"]):.2f}')
            bits = next((float(b[k]["bits_per_cell"]) for k in ("compact", "hashed") if k in b and "bits_per_cell" in b[k]), None)
            tiles.append(("Index memory (GB)", v, "; ".join(parts) + (f"; {bits:.1f} bits a value cell, not 64" if bits else "") + f'. {int(b["build"]["used_keys"]) / 2**30 * 100:.1f}% of the 4<sup>15</sup> keys are used.'))
        for w in WORKS[:2]:
            x0 = t(d, BASE, w, threads[0])
            cands = [(t(d, v, w, threads[0]), v) for v in PACKED if t(d, v, w, threads[0])]
            if x0 and cands:
                x1, v = max(cands, key=lambda c: c[0]["mqs"])
                tiles.append((f"{w[2].split(',')[0]}, 1 thread", f"{x1['mqs'] / x0['mqs']:.2f}× faster", f"{x0['mqs']:.2f} → {x1['mqs']:.2f} M seed lookups/s ({v[0]} keys)"))
        bad = mismatches(d)
        tiles.append(("Same answers", "yes" if not bad else "NO", "the hash of every emitted (position, reference, distance, ties) per query agrees across all variants" if not bad else "differ: " + ", ".join(bad)))
        a('<div class="tiles">' + "".join(f'<div class="tile"><div class="k">{k}</div><div class="v{" fail" if v == "NO" else ""}">{v}</div><div class="s">{s}</div></div>' for k, v, s in tiles) + "</div>")

        shown = [v for v in VARIANTS if v[4] and (d,) + v[:3] in loads]
        a('<h3>Lookups per second (prefetch on)</h3><div class="legend">' + "".join(f'<span><span class="sw" style="background:{v[4]}"></span>{html.escape(v[3])}</span>' for v in shown) + "</div>")
        groups = []
        for w in WORKS[:3]:
            for th in threads:
                vals = [(v, t(d, v[:3], w, th)) for v in shown]
                if all(x for _, x in vals):
                    groups.append((f"{w[3].replace('<br>', ' ')} · {th} thr", [(v[4], x["mqs"]) for v, x in vals]))
        if groups:
            a(bar_chart(groups))

        a('<h3>Every variant</h3><p class="muted">M queries/s; the ratio is against the first row in the same column. "pf": the batch prefetch (the old blob prefetches the key block 24 lookups ahead; the packed index runs its pipeline).</p>')
        cols = [(w, th, pf) for w in WORKS for th in threads for pf in (0, 1)]
        a('<div class="scroll"><table><tr><th>Variant</th><th class="num">Load s</th><th class="num">RSS GB</th>' + "".join(f'<th class="num">{w[3]}<br>{th}t{" pf" if pf else ""}</th>' for w, th, pf in cols) + "</tr>")
        for v in VARIANTS:
            key = (d,) + v[:3]
            if key not in loads:
                continue
            l = loads[key]
            a(f'<tr class="{"base" if v[:3] == BASE else ""}"><td>{html.escape(v[3])}</td><td class="num">{l["load_s"]:.1f}</td><td class="num">{l["rss_gb"]:.2f}</td>')
            for w, th, pf in cols:
                x, x0 = t(d, v[:3], w, th, pf), t(d, BASE, w, th, pf)
                if not x:
                    a('<td class="num">–</td>')
                    continue
                r = x["mqs"] / x0["mqs"] if x0 else 1
                rc = "up" if r > 1.05 else "down" if r < 0.95 else "muted"
                a(f'<td class="num">{x["mqs"]:.2f}' + ("" if v[:3] == BASE else f' <span class="{rc}">{r:.2f}×</span>') + "</td>")
            a("</tr>")
        a("</table></div>")
        x, xv = t(d, BASE, WORKS[0], threads[0], 0), t(d, BASE, WORKS[3], threads[0], 0)
        if x and xv:
            a(f'<p class="muted">Seeds emitted per lookup: {x["values"] / x["q"]:.1f}; values per key looked up: {xv["values"] / xv["q"]:.1f}.</p>')

    a(STRATEGIES)
    sp = os.path.join(rdir, "strategies_d2.tsv")
    if os.path.exists(sp):
        a('<div class="scroll"><table><tr><th>Packed compact, D2, 1 thread, prefetch on</th><th class="num">150 bp stream</th><th class="num">150 bp per read</th><th class="num">10 kb per read</th></tr>')
        cfg = {}
        for f in rows(sp):
            if f[2] == "pf=1":
                cfg.setdefault(f[0], {})[f[1]] = float(f[3])
        names = {"onepass=0 lines=8 step=8": "As first built (two-pass flank scan, 8 value lines, step 8)", "onepass=1 lines=8 step=8": "One-pass flank scan", "onepass=1 lines=16 step=8": "+ 16 value lines prefetched", "onepass=1 lines=32 step=8": "+ 32 value lines", "onepass=1 lines=32 step=4": "32 lines, step 4", "onepass=1 lines=32 step=12": "32 lines, step 12"}
        for c, m in cfg.items():
            a(f'<tr><td>{names.get(c, c)}</td>' + "".join(f'<td class="num">{m.get(s, 0):.2f}</td>' for s in ("reads", "perread", "long")) + "</tr>")
        a("</table></div>")
    a(SETUP)
    a("</div>")
    open(out, "w").write("\n".join(h))


def bar_chart(groups):
    """Horizontal bars, one scale; one bar per variant in each group, the ratio against the first."""
    w, left, right = 860, 190, 80
    bh, gap, gg = 13, 2, 16
    n = len(groups[0][1])
    hgt = len(groups) * (n * (bh + gap) + gg) + 10
    mx = max(v for _, bars in groups for _, v in bars) * 1.05
    sx = lambda v: left + (w - left - right) * v / mx
    s = [f'<div class="scroll"><svg viewBox="0 0 {w} {hgt}" width="100%" style="max-width:{w}px;min-width:560px" role="img" aria-label="Lookups per second by variant">']
    y = 6
    for label, bars in groups:
        s.append(f'<text x="{left - 10}" y="{y + (n * (bh + gap)) / 2 + 3}" text-anchor="end" font-size="12.5">{html.escape(label)}</text>')
        base = bars[0][1]
        for i, (col, v) in enumerate(bars):
            yy = y + i * (bh + gap)
            s.append(f'<rect x="{left}" y="{yy}" width="{max(sx(v) - left, 1):.1f}" height="{bh}" rx="3" fill="{col}"><title>{v:.2f} M queries/s</title></rect>')
            s.append(f'<text x="{sx(v) + 6:.1f}" y="{yy + bh - 2}" font-size="11.5" class="mono">{v:.2f}{"" if i == 0 else f"  {v / base:.2f}×"}</text>')
        y += n * (bh + gap) + gg
    s.append(f'<line class="axis" x1="{left}" x2="{left}" y1="0" y2="{y - gg + 4}"/></svg></div>')
    return "".join(s)


CSS = """
:root{--bg:#f6f7f4;--surface:#ffffff;--fg:#15181a;--muted:#5a6066;--line:#dcdfd9;--old:#eb6834;--new:#2a78d6;--hash:#1baf7a;--good:#0b7a35;--bad:#b4232c;
--sans:"IBM Plex Sans",system-ui,sans-serif;--mono:"IBM Plex Mono",ui-monospace,monospace;--display:"IBM Plex Sans Condensed","IBM Plex Sans",sans-serif}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#131514;--surface:#1b1e1c;--fg:#eceeea;--muted:#a3aaa4;--line:#313632;--old:#d95926;--new:#3987e5;--hash:#199e70;--good:#3fbf6a;--bad:#ef6b6b;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#131514;--surface:#1b1e1c;--fg:#eceeea;--muted:#a3aaa4;--line:#313632;--old:#d95926;--new:#3987e5;--hash:#199e70;--good:#3fbf6a;--bad:#ef6b6b;color-scheme:dark}
body{background:var(--bg);color:var(--fg);font:15px/1.55 var(--sans);margin:0}
.wrap{max-width:1100px;margin:0 auto;padding-block:28px 64px;padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px))}
h1{font:600 30px/1.15 var(--display);margin:0 0 6px;text-wrap:balance} h2{font:600 21px/1.25 var(--display);margin:40px 0 6px;text-wrap:balance} h3{font:600 15px/1.3 var(--sans);margin:22px 0 6px}
p{max-width:76ch;margin:6px 0} .lede{color:var(--muted);font-size:15.5px} .muted{color:var(--muted)} .mono,td.num,th.num{font-family:var(--mono);font-variant-numeric:tabular-nums}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:12px;margin:14px 0}
.tile{background:var(--surface);border:1px solid var(--line);border-radius:8px;padding:14px 16px;min-width:0}
.tile .k{font-size:12px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted)} .tile .v{font:600 27px/1.15 var(--display);margin:4px 0} .tile .s{font-size:13px;color:var(--muted)}
.scroll{overflow-x:auto;margin:8px 0 4px} table{border-collapse:collapse;font-size:13px;width:100%;background:var(--surface)}
th,td{padding:6px 9px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap} th{font-weight:600;font-size:11.5px;color:var(--muted);letter-spacing:.03em;vertical-align:bottom}
td.num,th.num{text-align:right} tr.base td{font-weight:600}
.up{color:var(--good)} .down{color:var(--bad)} .ok{color:var(--good);font-weight:600} .fail{color:var(--bad);font-weight:600}
.legend{display:flex;gap:16px;flex-wrap:wrap;font-size:13px;color:var(--muted);margin:4px 0}.sw{display:inline-block;width:11px;height:11px;border-radius:3px;margin-right:6px;vertical-align:-1px}
svg text{fill:var(--fg);font-family:var(--sans)} svg .axis{stroke:var(--line)} ul{max-width:80ch;padding-left:20px} li{margin:4px 0}
"""

STRATEGIES = """<h2>What was changed, and what was tried</h2><ul>
<li><b>Packed values.</b> A value takes ref_bits + pos_bits (17 + 24 = 41 bits at D1), not 64; flank cells stay 32-bit and byte-aligned so the flank scan reads plain integers; every value can be read on its own; the key table's offsets are unchanged.</li>
<li><b>Compact key table.</b> A presence bit per key (448 keys and a rank count per 64-byte line) and a 24-bit start per used key, instead of a 16-bit cell for every one of the 4<sup>15</sup> keys (2.68 GB). A hit reads two lines, one after the other.</li>
<li><b>Hashed key table</b> (added after the first results). The used keys in a table of 64-byte lines, 8 entries a line at load ≤ 0.75, so a hit reads one line; 10.7 bytes a used key. A bit per line marks that keys homed there overflowed, so most misses read one line (random-key lookups 17 → 27 M/s at D1 with it); misses still cost more than in the compact table. Built from the compact file in 0.5 s, sharing its values. Its line count is any number (a key's home line is key × lines / 4<sup>15</sup>), not a power of two: D3's 50.43 M used keys took 1.07 GB of lines as a power of two and take 0.54 GB now (same lookups, 3.61 → 3.51 M/s per read at the lower load).</li>
<li><b>Prefetch pipeline.</b> The key table 24 lookups ahead, a compact table's start 16 ahead, the values 8 ahead. The old blob prefetches the key block 24 ahead only, which never starts inside a 150 bp read (~24 k-mers); the packed pipeline still overlaps a read's own lookups.</li>
<li><b>Huge pages</b> (madvise) for the large arrays, also in the old blob's loader: 1.2–1.4× for the old blob alone. <span class="mono">FLEXMAP_THP=0</span> turns them off (the baseline rows).</li>
<li><b>Parallel load</b> of the packed file into a huge-page buffer; mmap is also offered.</li>
<li><b>Two FMKeys fixes:</b> the occurrence count saturates instead of wrapping at 65,536, and a key whose cells would pass its block's 16-bit offset is left out instead of corrupting the block.</li>
<li><b>Tried, no gain on D2</b> (table below): a one-pass flank scan, prefetching 16 or 32 value lines instead of 8, prefetch steps of 4 or 12, and a faster value decoder. D2 is bound by emitting its ~150 seeds a lookup (about 7 ns each in the consumer), which every variant pays alike.</li></ul>"""

SETUP = """<h2>Setup</h2><ul>
<li>Index as flexalign builds it: K 31, C 15, F 16, closed syncmer S 4, 16 keys per block, header threshold 2, keys over 2,048 values left out.</li>
<li>Queries are the true k-mers of contiguous sequence: reads drawn from the reference (random strand, 1% substitutions); every syncmer k-mer of a read, in read order, is a query, with its key and flank as the build makes them. <b>Per read</b>: one lookup call per read, as an aligner takes a read at a time. <b>Stream</b>: all reads' k-mers in one call per thread. <b>Random keys</b>: uniform over 2<sup>30</sup>.</li>
<li>The consumer hashes every emitted seed (position, reference, distance, ties) into a per-query checksum; the sums are compared across variants. The unit tests compare every key of a test index field by field.</li>
<li>Indexes on /dev/shm (RAM) or on SSD (D1 and D3 hashed), so load times are mostly memory speed. Each variant runs in its own process; threads take contiguous parts of the reads.</li>
<li>Machine: Intel Core i9-13900HX (AVX2, no AVX-512), 12 logical CPUs, 31 GB RAM shared with other work, Linux 6.17, rustc nightly 1.99. Some loads ran under memory pressure, so load times vary; runs that swapped during lookups were repeated.</li></ul>"""

if __name__ == "__main__":
    main()
