# SHiP hands-on: bash, git & ROOT

A short, hands-on tutorial for CERN summer students. It's a **guided walkthrough** of
the everyday analysis commands — **plotting**, **cutting**, and **fitting** — on a real
SHiP-like sample (a Heavy Neutral Lepton decaying to two daughters).

**Nothing to solve.** Every cell is already written: you run them and watch what
happens, then tweak a number here and there to see the effect. It all runs in a
**Jupyter notebook on SWAN**, with **git** to get the code and a little **bash** to
get around.

> **Every command you might need — bash, git, ROOT — is listed in [`TUTORIAL.md`](TUTORIAL.md).** Keep it open in a tab as you go.

---

## 1. Open SWAN (nothing to install)

1. One time only: visit [cernbox.cern.ch](https://cernbox.cern.ch) and log in (activates your storage).
2. Go to [swan.cern.ch](https://swan.cern.ch) → log in → **Start my Session** (default stack is fine).
3. Open the **terminal** (top-right icon).

## 2. Get the code (bash + git)

```bash
git clone https://github.com/vasilllis/ship-tutorial.git
cd ship-tutorial
ls -lh
```

`git clone` brings the project to SWAN; `ls` shows what's inside.

## 3. Run the notebook (ROOT)

The data lives on **EOS** at `/eos/experiment/ship/user/matclim/` — you have access, so
there's nothing to download (the data is **not** in this git repo).

In SWAN, open **`tutorial.ipynb`** and run the cells top to bottom (Shift+Enter). It
builds up the real analysis from `reference/` one small step at a time:

1. require two daughters (an RDataFrame `Filter`, with a cut-flow report),
2. plot the **displaced decay vertex**,
3. reconstruct the **invariant mass** (sum the daughter 4-momenta) — a peak at 1 GeV,
4. add detector **smearing** and watch the mass peak broaden,
5. reconstruct the **decay vertex** from the two tracks (the DOCA method — provided as a
   small toolbox you just call),
6. measure the **vertex resolution** vs smearing, with a built-in closure check at σ=0.

Every cell is provided; where you see **▸ Try**, change a value and re-run.

## 4. (optional) Save what you changed (git)

```bash
git add -A
git commit -m "Played with the tutorial"
```

---

## Files

```
ship-tutorial/
├── tutorial.ipynb     <- run this on SWAN
├── TUTORIAL.md        <- every command you might need (bash/git/ROOT)
├── reference/         <- the full analysis (vertex reco + resolution study), for the curious
├── README.md  SETUP.md  PUSH_TO_GITHUB.md
```

The full analysis (DOCA vertex reconstruction, five observables, the complete
resolution study) is in `reference/llp_simple_analysis.py` / `.C`, with a detailed
write-up in `reference/README.md`.
