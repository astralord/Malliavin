# Malliavin
Malliavin calculus on Haskell

Example: plot 3 trajectories of Brownian motions on the interval [0, 1] with discretization n = 500
```haskell
import Graphics.Plot

let lastT = 1
let n = 500
let b1 = brownianMotion 1 n lastT
let b2 = brownianMotion 2 n lastT
let b3 = brownianMotion 3 n lastT
mplot [linspace n (0, lastT), b1, b2, b3]
```
![alt tag](https://github.com/Quanteeks/Malliavin/blob/master/Brownian%20motions.png)
