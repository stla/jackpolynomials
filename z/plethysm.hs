-- plethysm' :: 
--   (Eq a, AlgField.C a) 
--   => Int 
--   -> [(Partition, RatioOfSprays a)] 
--   -> [(Partition, RatioOfSprays a)] 
--   -> ParametricSpray a
-- plethysm' n f g = 
--   foldl' 
--     (\spray (lambda, c) -> 
--       spray ^+^ c *^ productOfSprays (map plethysm_g lambda)) 
--     zeroSpray f
--   where 
--     plethysm_mu mu k = psPolynomial n (map (*k) mu)
--     q = lone 1
--     t = lone 2
--     plethysm_g k = 
--       foldl' (\spray (mu, r) -> spray ^+^ (changeVariables r [q, t^**^k]) *^ plethysm_mu mu k) zeroSpray g

-- plethysm :: 
--   (Eq a, AlgField.C a) 
--   => Int 
--   -> Map Partition (RatioOfSprays a)
--   -> (Int, Map Partition (Spray a -> RatioOfSprays a)) 
--   -> ParametricSpray a
-- plethysm n f (i, g) = 
--   DM.foldlWithKey 
--     (\spray lambda c -> 
--         spray ^+^ c *^ productOfSprays (map plethysm_g lambda)) 
--     zeroSpray f
--   where 
--     plethysm_mu mu k = psPolynomial n (map (*k) mu)
--     t = lone' i
--     plethysm_g k = 
--       DM.foldlWithKey 
--         (\spray mu d -> spray ^+^ (d (t k) *^ plethysm_mu mu k)) 
--           zeroSpray g

