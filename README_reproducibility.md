# Reproducibility Guide (Certificates + Build)

This package contains all code and tables to **reproduce the numeric certificates** and **compile the manuscript**.

## 1. Environment
- Python 3.10+
- `mpmath` (for high-precision arithmetic), `pandas` (for CSV handling)

Install:
```bash
python3 -m pip install mpmath pandas
```

## 2. Generate Certificates
### SU(2)
```bash
python3 routeB_certificates.py
```
Creates:
- `routeB_certificates.csv`
- `routeB_certificates_table.tex`

### SU(3)
```bash
python3 routeB_su3_certificates.py
```
Creates:
- `routeB_su3_certificates.csv`
- `routeB_su3_certificates_table.tex`

Acceptance criteria per row: `S(tau) < 1` and `kappa_A > 0`.

## 3. Compile the Manuscript
Use `latexmk` or `pdflatex` (run twice if needed):
```bash
latexmk -pdf main.tex
# or
pdflatex main.tex && pdflatex main.tex
```

Make sure the generated `*_table.tex` files are in the LaTeX working directory (they are included via `\input{...}`).

## 4. Notes
- The **benchmarks table** in the front-matter reports canonical constants (`C_H(a)`, `lambda_1`) and a representative SU(3) row from the CSV.
- Both Routes A (geometric) and B (constructive) are cross-referenced from the **Main Theorem** for immediate verification.