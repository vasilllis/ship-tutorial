# SHiP hands-on: bash, git & ROOT

A short, hands-on tutorial for CERN summer students. It's a **guided walkthrough** of
the everyday analysis commands — **plotting**, **cutting**, and **fitting** — on a real
SHiP-like sample (a Heavy Neutral Lepton decaying to two daughters).

**Nothing to solve.** Every cell is already written: you run them and watch what
happens, then tweak a number here and there to see the effect. It all runs in a
**Jupyter notebook on SWAN**, with **git** to get the code and a little **bash** to
get around.

---

## 1. Open SWAN (nothing to install)

1. One time only: visit <https://cernbox.cern.ch> and log in (activates your storage).
2. Go to <https://swan.cern.ch> → log in → **Start my Session** (default stack is fine).
3. Open the **terminal** (top-right icon).

## 2. Get the code (bash + git)

```bash
git clone https://github.com/$username$/ship-tutorial.git
cd ship-tutorial
ls -lh
```

`git clone` brings the project to SWAN; `ls` shows what's inside.

## 3. Run the notebook (ROOT)

The data lives on **EOS** at `/eos/experiment/ship/user/matclim/` — you have access, so
there's nothing to download (the data is **not** in this git repo).

In SWAN, open **`tutorial.ipynb`** and run the cells top to bottom (Shift+Enter). You'll:

- **plot** the decay vertex (the HNL decays tens of metres downstream) and the
  reconstructed **invariant mass** (a sharp peak at 1 GeV),
- **cut** on the fiducial decay volume and see how many events survive,
- **fit** a Gaussian to the mass peak to measure it,
- (optional) watch the peak **broaden** when detector resolution is added.

Where you see **▸ Try**, change a value and re-run — that's the exploring part.

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
├── reference/         <- the full analysis (vertex reco + resolution study), for the curious
├── README.md  SETUP.md  PUSH_TO_GITHUB.md
```

The full analysis (DOCA vertex reconstruction, five observables, the complete
resolution study) is in `reference/llp_simple_analysis.py` / `.C`, with a detailed
write-up in `reference/README.md`.
