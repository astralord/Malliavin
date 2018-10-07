module Malliavin where

import Data.Random.Distribution.Normal
import Numeric.LinearAlgebra.HMatrix as H
import Data.Ratio

data L2Map a = L2Map {eval :: a -> Double,
                      norm_l2 :: a -> Double}
l2_1 = L2Map {eval = (\t -> 1), norm_l2 = id}
l2id = L2Map {eval = id, norm_l2 = (\x -> x ^ 3 / 3)}

data L1Map a = L1Map {norm_l1 :: a -> Double}
l1_1 = L1Map {eval = (\t -> 1), norm_l1 = id}
l1id = L1Map {eval = id, norm_l1 = (\x -> x ^ 2 / 2)}


-- | 'hermite' is an infinite list of Hermite polynomials for given x
hermite :: (Enum a, Num a) => a -> [a]
hermite x = s
  where s@(_:ts) = 1 : x : zipWith3 (\hn2 hn1 n1 -> x * hn1 - n1 * hn2) s ts [1..]


-- | 'gaussianProcess' samples from Gaussian process
gaussianProcess :: Seed                   -- random state
                -> Int                    -- number of samples m
                -> Int                    -- number of dimensions n
                -> ((Int, Int) -> Double) -- function that maps indices of the matrix into dot products of its elements
                -> [Vector Double]        -- m n-th dimensional samples of Gaussian process
gaussianProcess seed 0 _ _ = []
gaussianProcess seed m 0 _ = replicate m $ vector []
gaussianProcess seed m n dotProducts = toRows $ gaussianSample seed m mean cov_matrix
  where mean = n |> repeat 0
        cov_matrix = H.sym $ (n><n) $ map (\i -> dotProducts (quot i n, rem i n)) [0..]


-- | 'second' function takes a list and gives each second element of it
second (x:y:xs) = y : second xs
second _ = []


-- | 'coinTossExpansion' is a Wiener chaos expansion for coin-toss r.v. to n-th element
coinTossExpansion :: Int    -- number of elements in the sum
                  -> Double -- gaussian random variable
                  -> Double -- the sum
coinTossExpansion n xi = sum (take n $ 0.5 : zipWith (*) fn (second $ hermite xi))
  where fn = 1.0 / (sqrt $ 2 * pi) : zipWith (\fn1 k -> -fn1 * k / (k ^ 2 + 3 * k + 2)) fn [1, 3..]


-- | 'coinTossSequence' is a coin-toss sequence of given size
coinTossSequence :: Seed  -- random state ω
                 -> Int   -- size of resulting sequence m
                 -> [Int] -- coin-toss sequence of size m
coinTossSequence seed m = map (round.coinTossExpansion 100) (toList nvec)
  where nvec = gaussianProcess seed m 1 (\(i,j) -> 1) !! 0


-- | 'mapOverInterval' map function f over the interval [0, T]
mapOverInterval :: (Fractional a) => Int      -- n, size of the output list
                                  -> a        -- T, end of the time interval
                                  -> (a -> b) -- f, function that maps from fractional numbers to some abstract space
                                  -> [b]      -- list of values f(t), t \in [0, T]
mapOverInterval n endT fn = [fn $ (endT * fromIntegral i) / fromIntegral (n - 1) | i <- [0..(n-1)]]


type ItoIntegral = Seed         -- ω, random state
                -> Int          -- n, sample size
                -> Double       -- T, end of the time interval
                -> L2Map Double -- h, L2-function
                -> [Double]     -- list of values sampled from Ito integral


-- | 'itoIntegral'' trajectory of Ito integral on the time interval [0, T]
itoIntegral' :: ItoIntegral
itoIntegral' seed n endT h = 0 : (toList $ gp !! 0)
  where gp = gaussianProcess seed 1 (n-1) (\(i, j) -> norm_l2 h $ fromIntegral (min i j + 1) * t)
        t = endT / fromIntegral n


-- | 'itoIntegral' faster implementation of itoIntegral' function
itoIntegral :: ItoIntegral
itoIntegral seed 0 _    _ = []
itoIntegral seed n endT h = scanl (+) 0 increments
  where increments = toList $ sigmas * gaussianVector
        gaussianVector = flatten $ gaussianSample seed (n-1) (vector [0]) (H.sym $ matrix 1 [1])
        sigmas = fromList $ zipWith (\x y -> sqrt(x-y)) (tail hnorms) hnorms
        hnorms = mapOverInterval n endT $ norm_l2 h


-- | 'brownianMotion' trajectory of Brownian motion a.k.a. Wiener process on the time interval [0, T]
brownianMotion :: Seed     -- ω, random state
               -> Int      -- n, sample size
               -> Double   -- T, end of the time interval
               -> [Double] -- list of values sampled from Brownian motion
brownianMotion seed n endT = itoIntegral seed n endT l2_1


type DiffusionProcess = Seed         -- ω, random state
                     -> Int          -- n, sample size
                     -> Double       -- T, end of the time interval
                     -> L1Map Double -- function \mu with L1 norm
                     -> L2Map Double -- function \sigma with L2 norm
                     -> [Double]     -- list of values sampled from stochastic diffusion process X(t)

-- | 'diffusionProcess' stochastic process, defined by SDE: dX(t) = \mu(t)dt + \sigma(t)dB(t)
diffusionProcess :: DiffusionProcess
diffusionProcess seed n endT mu sigma = zipWith (+) trend noise
  where trend = mapOverInterval n endT $ norm_l1 mu
        noise = itoIntegral seed n endT sigma


-- | 'multipleItoIntegral' multiple m-th Ito integral of m-th tensor power of h
multipleItoIntegral :: Seed         -- ω, random state
                    -> Int          -- m, tensor power
                    -> Int          -- n, sample size
                    -> Double       -- T, end of the time interval
                    -> L2Map Double -- function h with L2 norm
                    -> [Double]     -- list of values sampled from multiple m-th Wiener-Ito integral of h^(⊗m)
multipleItoIntegral seed m n endT h = zipWith fn itoTrajectory hnorm
  where itoTrajectory = itoIntegral seed n endT h
        hnorm = mapOverInterval n endT $ norm_l2 h
        fn :: Double -> Double -> Double
        fn _ 0 = 0
        fn x y = (y ^ m) * (hermite (x / y) !! m)

multipleBrownianMotion :: Seed     -- ω, random state
               -> Int      -- m, tensor power
               -> Int      -- n, sample size
               -> Double   -- T, end of the time interval
               -> [Double] -- list of values sampled from Brownian motion
multipleBrownianMotion seed m n endT = multipleItoIntegral seed m n endT l2_1


-- | 'totalVariation' calculates total variation of the sample process X(t)
totalVariation :: [Double] -- process trajectory X(t)
               -> Double   -- V(X(t))
totalVariation x@(_:tx) = sum $ zipWith (\x xm1 -> abs $ x - xm1) tx x


-- | 'quadraticVariation' calculates quadratic variation of the sample process X(t)
quadraticVariation :: [Double] -- process trajectory X(t)
                   -> Double   -- [X(t)]
quadraticVariation x@(_:tx) = sum $ zipWith (\x xm1 -> (x - xm1) ^ 2) tx x
