# [**SVDFR**](https://gokulhari.github.io/svdfr/) 

## Singular value decomposition of the frequency response operator using [Chebfun](http://www.chebfun.org).

[**SVDFR**](https://gokulhari.github.io/svdfr/) can solve for singular values of frequency responses of a subset of linear partial differential equations in time and one other spatial variable. The principal singular value computed in this manner is the largest possible amplification of the output for a unit input. The associated functions of the principal singular value give the shape of the input that can lead to the optimal output.

[Lieu and Jovanovic](https://doi.org/10.1016/j.jcp.2013.05.010) originally started this endeavor in using [Chebfun's](http://www.chebfun.org) automatic collocation in frequency response calculations. They used an integral reformulation to derive a well-conditioned numerical scheme to compute frequency responses. Over the years, [Chebfun](http://www.chebfun.org) has evolved substantially, and the original version of SVDFR is no more maintained. In contrast, the version on this site uses [Chebfun's](http://www.chebfun.org) built-in well-conditioned [ultraspherical method](https://doi.org/10.1137/120865458), and a robust feedback interconnected system (see the [accompanying paper](https://arxiv.org/abs/2005.04493)). Furthermore, the current version is updated to comply with the most recent release of [Chebfun](http://www.chebfun.org) (v5.7.0) and is designed to work directly with the [Chebfun](http://www.chebfun.org) syntax: *The end result as you will see in the examples is short and intuitive codes to compute frequency responses of relatively complex systems*.  

We also provide a routine to compute the $H_\infty$ norm of a linear system (the largest energy amplification among inputs of all temporal frequencies). Most users familiar with [Chebfun](http://www.chebfun.org) can start using [**SVDFR**](https://gokulhari.github.io/svdfr/) effortlessly. A more detailed description with examples can be found [here.](https://gokulhari.github.io/svdfr/)


 
