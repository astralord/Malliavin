import Data.Random.Distribution.Normal
import Numeric.LinearAlgebra.HMatrix as H

hermite :: (Enum a, Num a) => a -> [a]
hermite x = s
    where s@(_:ts) = 1 : x : zipWith3 (\hn2 hn1 n1 -> x * hn1 - n1 * hn2) s ts [1..]
 
gaussianProcess :: Seed -> Int -> Int -> ((Int, Int) -> Double) -> Matrix Double
gaussianProcess seed nSamples dim dotProducts = gaussianSample seed nSamples mean cov_matrix
    where mean = vector (replicate dim 0)
          cov_matrix = H.sym $ matrix dim (map (\i -> dotProducts (i `quot` dim, i `rem` dim)) $ take (dim * dim) [0, 1..])

second (x:y:xs) = y : second xs
second _ = []

coinTossExpansion :: Int -> Double -> Double
coinTossExpansion n xi = sum (take n $ 0.5 : zipWith (*) fn (second $ hermite xi))
    where fn = 1.0 / (sqrt $ 2 * pi) : zipWith ( \fn1 k -> -fn1 * k / ((k + 1) * (k + 2)) ) fn [1, 3..]


coinTossSequence :: Seed -> Int -> [Int]
coinTossSequence seed n = map (round.coinTossExpansion 100) (toList nvec)
    where nvec = (toColumns $ gaussianProcess seed n 1 (\(i,j)->1) ) !! 0
