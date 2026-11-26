#!/usr/bin/env bash
# Rebuild the proteomics panels that depend on the GO GSEA tables, from
# enrichment_within_compartment_20261125 (fgsea nPermSimple = 100000).
# Order: 02 (D), 04 (F, Supp F), 06 (H, I first pass), 12, 14, 15 (style
# matching), 16b (G), 17 (H, I final). One R process at a time; each step waits
# while 4 or more R processes are running on the machine. Logs in ../logs.
# Refuses an output root that holds anything besides scripts/ unless
# FIGURE_REBUILD_FORCE=true (the R scripts refuse to replace a file likewise).

set -u
RS="C:/Users/topohl/AppData/Local/Programs/R/R-4.5.1/bin/x64/Rscript.exe"
ROOT="S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/figure_rebuild_within_compartment_20261125"
SCRIPTS="$ROOT/scripts"
LOGS="$ROOT/logs"

other="$(ls -A "$ROOT" 2>/dev/null | grep -v '^scripts$')"
if [ -n "$other" ] && [ "${FIGURE_REBUILD_FORCE:-false}" != "true" ]; then
  echo "Output root is not empty: $ROOT (set FIGURE_REBUILD_FORCE=true to reuse it)." >&2
  exit 1
fi
mkdir -p "$LOGS"

r_count() {
  powershell.exe -NoProfile -Command "(Get-CimInstance Win32_Process -Filter \"Name='Rterm.exe' OR Name='Rscript.exe'\" | Measure-Object).Count" | tr -d '\r '
}
wait_for_r_slot() {
  while [ "$(r_count)" -ge 4 ]; do sleep 30; done
}

for s in 02_panelD_volcano 04_panelF_suppF_gsea_dotplots 06_panelH_scatter_panelI_heatmap \
         12_style_match_volcano 14_style_match_dotplots 15_style_match_scatter_heatmap \
         16b_final_Fig3G_gsea_lollipop 17_final_Fig3H_Fig3I; do
  wait_for_r_slot
  echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) $s"
  "$RS" "$SCRIPTS/$s.R" > "$LOGS/$s.log" 2>&1
  code=$?
  if [ "$code" -ne 0 ]; then echo "Step failed: $s (see $LOGS/$s.log)" >&2; exit "$code"; fi
done
echo "Rebuild completed: $ROOT"
