import Data.Random.Distribution.Normal
import Numeric.LinearAlgebra.HMatrix as H
import Graphics.Plot

-- | 'hermite' function returns infinite list of Hermite polynomials for given x
hermite :: (Enum a, Num a) => a -> [a]
hermite x = s
    where s@(_:ts) = 1 : x : zipWith3 (\hn2 hn1 n1 -> x * hn1 - n1 * hn2) s ts [1..]
 

-- | 'gaussianProcess' function takes random seed, number of samples, number of dimensions n and
-- function that takes indices and returns dot products of elements of the matrix with corresponding indices,
-- and returns Gaussian process as n-th dimensional vector
gaussianProcess :: Seed -> Int -> Int -> ((Int, Int) -> Double) -> Vector Double
gaussianProcess seed nSamples dim dotProducts = flatten $ gaussianSample seed nSamples mean cov_matrix
    where mean = vector (replicate dim 0)
          cov_matrix = H.sym $ matrix dim (map (\i -> dotProducts (i `quot` dim, i `rem` dim)) $ take (dim * dim) [0, 1..])


-- | 'second' function takes a list and returns each second element of it
second (x:y:xs) = y : second xs
second _ = []


-- | 'coinTossExpansion' takes an integer parameter n and gaussian r.v. xi,
-- and returns its Wiener chaos expansion to n-th element with parameter xi
coinTossExpansion :: Int -> Double -> Double
coinTossExpansion n xi = sum (take n $ 0.5 : zipWith (*) fn (second $ hermite xi))
    where fn = 1.0 / (sqrt $ 2 * pi) : zipWith ( \fn1 k -> -fn1 * k / ((k + 1) * (k + 2)) ) fn [1, 3..]


-- | 'coinTossSequence' takes random seed and integer n and returns a coin-toss sequence of size n
coinTossSequence :: Seed -> Int -> [Int]
coinTossSequence seed n = map (round.coinTossExpansion 100) (toList nvec)
    where nvec = gaussianProcess seed n 1 (\(i,j) -> 1)


itoIntegral :: Seed -> Int -> Double -> (Double -> Double) -> Vector Double
itoIntegral seed n lastT hnorm = mappend 0 $ gaussianProcess seed 1 (n - 1) (\(i,j) -> hnorm $ fromIntegral (min i j + 1) * t)
    where t = lastT / fromIntegral n

testFun :: Double -> Double
testFun t = 0.5 * t - sin (2 * pi * t) / (4 * pi)


-- | 'brownianMotion' function takes random seed, parameters n and t and
-- returns n-th dimensional vector, filled with B(s) for s = [0, ... , t],
-- where B(s) stands for Brownian motion (Wiener process)
brownianMotion :: Seed -> Int -> Double -> Vector Double
brownianMotion seed n lastT = itoIntegral seed n lastT id

