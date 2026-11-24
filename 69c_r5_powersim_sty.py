#!/usr/bin/env python
"""R5 D2a (PREREG_r5.md): the frozen R2 simulator (63a3) with three weaker scenarios added. Same job syntax as 63a3:
CLADE HOG SCEN REP OUTDIR; OUTDIR/../_fit links to the R2 calibration, so no recalibration happens."""
import os, sys, importlib.util
B = os.path.expanduser("~/deepsea_conv")
sp = importlib.util.spec_from_file_location("ps", f"{B}/scripts/63a3_powersim_sty.py"); ps = importlib.util.module_from_spec(sp); sp.loader.exec_module(ps)
ps.SCEN.update(K85=("relax", 0.85, 0, 0, 0), K90=("relax", 0.90, 0, 0, 0), K75h=("relax", 0.75, 1, 0, 0))
if __name__ == "__main__":
    if not os.path.exists(os.path.join(os.path.dirname(sys.argv[5].rstrip("/")), "_fit", "scale.txt")): sys.exit("no R2 calibration")
    ps.main()
