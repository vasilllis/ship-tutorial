# SHiP hands-on: bash, git & ROOT

A short, hands-on tutorial for CERN summer students. You'll reconstruct a
**Heavy Neutral Lepton (HNL)** — a long-lived particle that decays to two
daughters at a displaced vertex — and study how detector resolution affects the
measurement. Everything runs in a **Jupyter notebook on SWAN**, with **git** to
keep your work and a little **bash** to get around.

It's deliberately simple: one notebook, the real analysis data, two `TODO` cells.

---

## 1. Open SWAN (nothing to install)

1. One time only: visit <https://cernbox.cern.ch> and log in (activates your storage).
2. Go to <https://swan.cern.ch> → log in → **Start my Session** (default stack is fine).
3. Open the **terminal** (top-right icon).

## 2. Get the repo (bash + git)

```bash
git clone https://github.com/<ORG>/ship-tutorial.git
cd ship-tutorial

ls -lh            # look around
ls -lh /eos/user/m/matclim/SHiPsim/SS_2026_data/   # the HNL sample lives on EOS
```

`clone` brought the project to SWAN; `ls` shows what's inside.

## 3. Do the analysis (ROOT, in the notebook)

> The notebook reads the sample directly from EOS
> (`/eos/user/m/matclim/SHiPsim/SS_2026_data`). You have EOS access, so there's
> nothing to download — the data is **not** stored in this git repo.

In SWAN, open **`tutorial.ipynb`** and run the cells top to bottom (Shift+Enter).
You fill in **two TODO cells**:

- **TODO 1** — reconstruct the transverse momentum `pt`.
- **TODO 2** — print the mass resolution (StdDev) at each smearing level.

Along the way you'll reconstruct the HNL mass (a sharp peak at **1 GeV**), see the
**displaced decay vertex** at tens of metres, and watch the mass peak **broaden**
as detector resolution worsens.

## 4. Save your work (git)

```bash
git add tutorial.ipynb
git commit -m "Complete the TODOs"
git log --oneline
```

---

## Stretch goals

- **Branch & merge.** Try a different smearing range on a branch, then merge it.
- **Make (and resolve) a conflict** on the same line — the git skill that matters most.
- **Go deeper.** The full analysis (vertex reconstruction via DOCA, vertex-resolution
  study, five observables) is in `reference/` — `llp_simple_analysis.py` / `.C`, with a
  detailed write-up in `reference/README.md`.

## Files

```
ship-tutorial/
├── tutorial.ipynb     <- start here
├── solutions.ipynb    <- the completed version
└── reference/         <- full original analysis (.C, .py, write-up)
└── reference/         <- the full original analysis (optional, for the curious)
```

The dataset is a simulated HNL of mass 1 GeV decaying to two daughters. Have fun.
