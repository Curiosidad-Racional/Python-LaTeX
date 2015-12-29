# Python-LaTeX

Python-LaTeX is a program written in python that process .ptx
files with latex code and, between this, python code
embedded with '««' '»»' delimiters. Python-LaTeX output
latex file.

## Dependences

- [TeX Live 2013](http://www.tug.org/texlive/) - `$ sudo apt-get install texlive-full`
- [Python3](https://www.python.org/downloads/) - `$ sudo apt-get install python3`
- [Sympy](http://www.sympy.org/en/download.html)
```
$ git clone git://github.com/sympy/sympy.git
$ cd sympy
$ sudo python3 setup.py install
```
- [SciPy](http://www.scipy.org/install.html) - `$ sudo apt-get install python3-scipy`
- [NumPy](http://www.scipy.org/scipylib/download.html) - `$ sudo apt-get install python3-numpy`

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


## Examples

### Basics

`examples/basics` folder.

Command line:

```
$ py2pdf -s Natural\ numbers.ptx 
% Python-Sympy -> LaTeX success.
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
```

View results:

```
$ evince Natural\ numbers.pdf
$ evince Natural\ numbers\ -\ Sol.pdf
```


### Plots

`examples/plots` folder.

Command line:

```
$ py2pdf -s Systems\ of\ linear\ equations.ptx 
% Python-Sympy -> LaTeX success.
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
```

View results:

```
$ evince Systems\ of\ linear\ equations.pdf
$ evince Systems\ of\ linear\ equations\ -\ Sol.pdf
```


### Cirtuits

`examples/circuits` folder.

Command line:

```
$ py2pdf -s Y-Delta\ circuit.ptx 
Ciclos independientes teóricos: 3 reales: 3
% Python-Sympy -> LaTeX success.
0m1s real, 0m1s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m1s real, 0m1s user, 0m0s sys
0m1s real, 0m1s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
0m0s real, 0m0s user, 0m0s sys
```

View results:

```
$ evince Y-Delta\ circuit.pdf
$ evince Y-Delta\ circuit\ -\ Sol.pdf
```

### Animations

`examples/animations` folder.

Command line:

```
$ py2pdf Two\ Charged\ Pendulums\ 3D.ptx 
Animation: Comenzando cálculos...
Cálculo explicito partículas: 0.
Cálculo Numérico partículas: [1, 2] 20m 12s .
Animation: Cálculos finalizados.
Procesando fotograma 301 de 301
Animation: Compilando código Latex del entorno picture 301...
latex procesó la página 301
dvips procesó la página 301
Animation: Compilado.
% Python-Sympy -> LaTeX success.
8m46s real, 8m44s user, 0m2s sys
0m1s real, 0m1s user, 0m0s sys
0m1s real, 0m1s user, 0m0s sys
0m6s real, 0m5s user, 0m0s sys
11m11s real, 10m51s user, 0m18s sys
```

View results:

```
$ acroread Two\ Charged\ Pendulums\ 3D.pdf
```

You need `acroread` to view animations

Install `acroread` on Ubuntu:

```
sudo add-apt-repository "deb http://archive.canonical.com/ raring partner"
sudo apt-get update
sudo apt-get install acroread
```