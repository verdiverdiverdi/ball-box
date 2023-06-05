# ball-box
volume of intersection for high dimensional centred balls and cubes

What is the volume of intersection of $B_n(1) \cap [-1/\sqrt{s}, 1/\sqrt{s}]^n$?
When $s \leq 1$ or $s \geq n$ it is the volume of the ball or the cube, respectively.
It was a fun adventure to find out when $s \in (1, n)$.
The question was originally answered in [1], and a survey given in [2] (which initially contained small errors).
I later learnt that [3] had derived a more general formula for non centred oblongs, rather than centred cubes.
The most recent revision of [2] is now correct.
An asymptotic estimate is also provided in [3], which requires $s > n/3$ for centred cubes.

I have chosen to use the method of Constales from [1], which represents the volume as an infinite sum of Fresnel integrals, using the [mpmath](https://mpmath.org/) package for its excellent selection of mathematical functions, and the ease with which it allows one to increase precision.
As such, your python will need access to it.

Given that I truncate an infinite sum, and use finite precision for Fresnel integrals, the volumes my script outputs are approximations.
For the case of $s > n/3$ one can use the asymptotic estimate as a ground truth (for large enough $n$, whatever that means) to check whether the number of terms and the precision I choose are sufficient.
For the case $s \in (1, 2)$, where only the spherical caps peak outside the cube, one can use the formulas of [4] to determine the accuracy of the approximation as $n$ increases.
Using an experimentally observed heuristic (via Monte--Carlo sampling methods in dimensions up to $n = 80$, not reproduced here) that the approximation becomes more accurate as $s$ grows, the $s \in (1, 2)$ case represents an experimental worst case, and I therefore assume that for $s \approx 1$ these approximations are good for $n \leq 220$, and for $s \geq 2$ these approximations are good for $n \leq 250$ (see `findAccurateRange` in `accuracyTests.py`).

Precomputed, high precision, values for the Fresnel integrals in a given dimension will be saved in a directory called `precomps` -- each will be about 3MB given current settings, so something to consider if planning on increasing the precision or number of terms, and approximating in many dimensions.

Future directions:
- More satisfying understanding of the required number of terms and precision for a given pair $(n, s)$ to achieve a sufficiently precise approximation.
- Implement the more general case given in [3].

To use these scripts, I open sage and then

```
attach("fresnel.py")
fresnelConstalesLog(n, s)
```

which returns $\textnormal{logvol}(B_n(1) \cap [-1/\sqrt{s}, 1/\sqrt{s}]^n)$ (in the natural log, of course).

[1]: Cecil C. Rousseau and Otto G. Ruehr. Problems and solutions. SIAM Review, 39(4):761–789, 1997.

[2]: Peter J. Forrester. Comment on “sum of squares of uniform random variables” by I. Weissman. https://arxiv.org/abs/1804.07861, 2018.

[3]: Yoshinori Aono and Phong Q. Nguyen. Random sampling revisited: Lattice enumeration with discrete pruning. In Jean-S\'{e}bastien Coron and Jesper Buus Nielsen, editors, EUROCRYPT 2017, Part II, volume 10211 of LNCS, pages 65–102. Springer, Heidelberg, April / May 2017.

[4]: S. Li , 2011. Concise Formulas for the Area and Volume of a Hyperspherical Cap. Asian Journal of Mathematics & Statistics, 4: 66-70. 10.3923/ajms.2011.66.70
