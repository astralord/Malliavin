hermite :: (Enum n, Num n) => n -> [n]
hermite x = s
    where s@(_:ts) = 1 : x : zipWith3 (\hn2 hn1 n1 -> x*hn1 - n1*hn2) s ts [1..]
