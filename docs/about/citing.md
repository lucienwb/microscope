# Citing

If microscope made a figure in your paper, or you used it to look at the
results that went into one, a citation helps other people find it:

> Wang, B. *microscope: a molecular structure and spectroscopy viewer for
> quantum chemistry*, version 0.1.0 (2026).
> <https://github.com/lucienwb/microscope>

```bibtex
@software{microscope,
  author  = {Wang, Bin},
  title   = {microscope: a molecular structure and spectroscopy viewer for
             quantum chemistry},
  version = {0.1.0},
  year    = {2026},
  url     = {https://github.com/lucienwb/microscope},
  license = {MIT}
}
```

The same information is in
[`CITATION.cff`](https://github.com/lucienwb/microscope/blob/main/CITATION.cff)
at the top of the repository, which is what GitHub's **Cite this repository**
button reads: it gives the citation in APA and BibTeX, ready to copy.

Please cite the version you used, which `scope --version` prints.

microscope only draws what your calculations produced. The programs that
produced them, and the methods and basis sets, have their own citations.
