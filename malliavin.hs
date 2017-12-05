import Data.Random
import Data.Random.Distribution.Normal
import Data.Matrix
import Control.Monad
import Data.Random.Distribution.MultivariateNormal

hermite :: (Enum a, Num a) => a -> [a]
hermite x = s
    where s@(_:ts) = 1 : x : zipWith3 (\hn2 hn1 n1 -> x * hn1 - n1 * hn2) s ts [1..]

second (x:y:xs) = y : second xs;
second _ = []

coin_toss_expansion :: Int -> Double -> Double
coin_toss_expansion n xi = sum (take n $ 0.5 : zipWith (*) fn (second $ hermite xi))
    where fn = 1.0 / (sqrt $ 2 * pi) : zipWith ( \fn1 k -> -fn1 * k / ((k + 1) * (k + 2)) ) fn [1, 3..]

rounded_coin_toss_expansion :: Int -> Double -> Int
rounded_coin_toss_expansion n xi = round $ coin_toss_expansion n xi

coin_toss_sequence :: Int -> RVar [Int]
coin_toss_sequence n = replicateM n $ liftM (rounded_coin_toss_expansion 100) (stdNormal)

main = sample $ coin_toss_sequence 20
