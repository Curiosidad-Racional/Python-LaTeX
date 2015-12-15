# Python-LaTeX

Python-LaTeX is a program written in python that process .ptx
files with latex code and, between this, python code
embedded with '««' '»»' delimiters. Python-LaTeX output
latex file.

## Dependences

- TeX Live 2013
- python3
- scipy
- numpy

## Setup

Add two lines to `.bashrc` file.

1. Declare system variable.
 `export PYTHON_LATEX=[/path/to/Python-LaTex]`
2. Add to PATH variable.
 `export PATH="$PATH:$PYTHON_LATEX/bin"`

## Basic usage

- .ptx file to .ltx file.
 `py2tex [.ptx file]`
- .ptx file to .pdf file.
 `py2pdf [.ptx file]`
- .ptx file to .pdf files with and without solutions.
 `py2pdf -s [.ptx file]`

## Usage without scripts (for example on Windows)

- .ptx file to .ltx file.
 `python3 /path/to/Python-Latex/python/pytex.py [.ptx file]`

## Basic example

.ptx file: `sample.ptx`

```
\begin{equation}
7! = ««factorial(7)»»
\end{equation}
```

Command line

```
$ python3 /path/to/Python-Latex/python/pytex.py sample.ptx 
% Python-Sympy -> LaTeX success.
```

Output file: `sample.ltx`

```
[... Latex preamble ...]

\begin{document}

\begin{equation}
7! = 5040
\end{equation}

\end{document}
```