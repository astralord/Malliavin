# Malliavin
Malliavin calculus on Haskell

Example: plot 3 trajectories of Brownian motions on the interval [0, 1] with discretization n = 500
```haskell
import Graphics.Plot

let endT = 1
let n = 500
let b1 = brownianMotion 1 n endT
let b2 = brownianMotion 2 n endT
let b3 = brownianMotion 3 n endT
mplot [linspace n (0, endT), b1, b2, b3]
```
![alt tag](https://github.com/Quanteeks/Malliavin/blob/master/Brownian%20motions.png)
