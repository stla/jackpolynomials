-- _hallInnerProductFromPScombinations :: 
--   forall b. (Eq b, AlgRing.C b)
--   => (Spray b -> Spray b -> Spray b)
--   -> Map Partition (Spray b)
--   -> Map Partition (Spray b)
--   -> Spray b         -- ^ parameter
--   -> Spray b 
-- _hallInnerProductFromPScombinations multabFunc psCombo1 psCombo2 alpha = 
--   AlgAdd.sum $ DM.elems
--     (merge dropMissing dropMissing (zipWithMatched f) psCombo1 psCombo2)
--   where
--     zlambda' :: Partition -> Spray b
--     zlambda' lambda = fromIntegral (zlambda lambda) 
--       AlgRing.* alpha AlgRing.^ (toInteger $ length lambda)
--     f :: Partition -> Spray b -> Spray b -> Spray b
--     f lambda coeff1 coeff2 = 
--       multabFunc (zlambda' lambda) (coeff1 AlgRing.* coeff2)

-- _hallInnerProductFromPScombinations :: 
--   forall b. (Eq b, AlgRing.C b)
--   => Map Partition b
--   -> Map Partition b
--   -> Spray b 
-- _hallInnerProductFromPScombinations psCombo1 psCombo2 = 
--   AlgAdd.sum $ DM.elems
--     (merge dropMissing dropMissing (zipWithMatched f) psCombo1 psCombo2)
--   where
--     alpha = lone 1 :: Spray b
--     zlambda' :: Partition -> Spray b
--     zlambda' lambda = zlambda lambda .^ (alpha ^**^ (length lambda))
--     f :: Partition -> b -> b -> Spray b
--     f lambda coeff1 coeff2 = 
--       (coeff1 AlgRing.* coeff2) *^ (zlambda' lambda) 

-- -- symbolicHallInnerProductFromPScombinations :: 
-- --   (Eq b, AlgField.C b) => Map Partition b -> Map Partition b -> Spray b
-- -- symbolicHallInnerProductFromPScombinations =
-- --   _symbolicHallInnerProductFromPScombinations 
-- --     (
-- --       _hallInnerProductFromPScombinations (^*^)
-- --     ) 
-- -- _symbolicHallInnerProductFromPScombinations :: 
-- --   (Eq b, AlgRing.C b) 
-- --   => (Map Partition (Spray b) -> Map Partition (Spray b) -> Spray b -> Spray b) 
-- --   -> Map Partition b -> Map Partition b -> Spray b
-- -- _symbolicHallInnerProductFromPScombinations func psCombo1 psCombo2 = func psCombo1' psCombo2' (lone 1)
-- --   where
-- --     psCombo1' = DM.map constantSpray psCombo1
-- --     psCombo2' = DM.map constantSpray psCombo2

-- jackInPSbasis ::
--   forall a. (Eq a, AlgField.C a) => Char -> Partition -> Map Partition (RatioOfSprays a)
-- jackInPSbasis which mu = 
--   DM.filter (/= AlgAdd.zero) 
--     (unionsWith (^+^) (DM.elems $ DM.mapWithKey combo_to_map jackCombo))
--   where
--     jackCombo = jackInMSPbasis which mu :: Map Partition (RatioOfSprays a)
--     combo_to_map lambda rOS = 
--       DM.map 
--         (\r -> (fromRational r :: a) AlgMod.*> rOS) (mspInPSbasis lambda)
-- --          (psCombination (msPolynomial (sum lambda) lambda))

-- skewJackPolTerm :: 
--   forall a. (Eq a, AlgField.C a) => Int -> Char -> Partition -> Partition -> Partition -> ParametricSpray a
-- skewJackPolTerm n which lambda mu nu = 
--   coeff *^ jackSymbolicPol n nu 'Q'
--   where
--     psCombo_lambda = jackInPSbasis 'Q' lambda
--     psCombo_mu = jackInPSbasis 'P' mu
--     psCombo_nu = jackInPSbasis 'P' nu
--     psCombo_mu_nu = psCombinationsProduct psCombo_mu psCombo_nu
--     coeff = 
--       evaluate (_hallInnerProductFromPScombinations psCombo_lambda psCombo_mu_nu) [lone 1 %//% unitSpray]  
--     -- psCombo_lambda = jackInPSbasis which lambda
--     -- psCombo_mu = jackInPSbasis which mu
--     -- psCombo_nu = jackInPSbasis which nu
--     -- psCombo_mu_nu = psCombinationsProduct psCombo_mu psCombo_nu
--     -- coeff = 
--     --   evaluate (_hallInnerProductFromPScombinations psCombo_lambda psCombo_mu_nu) [lone 1 %//% unitSpray] AlgField./ 
--     --     evaluate (_hallInnerProductFromPScombinations psCombo_nu psCombo_nu) [lone 1 %//% unitSpray]
--     -- hip = _hallInnerProductFromPScombinations (AlgRing.*) :: 
--     --   Map Partition (ParametricSpray a) -> Map Partition (ParametricSpray a) -> ParametricSpray a -> ParametricSpray a
--     -- alpha = lone 1 :: ParametricSpray a
--     -- coeff = hip psCombo_lambda psCombo_mu_nu alpha --AlgField./ hip psCombo_nu psCombo_nu alpha

-- skewJackPol :: 
--   (Eq a, AlgField.C a) => Int -> Char -> Partition -> Partition -> ParametricSpray a
-- skewJackPol n which lambda mu = 
--   sumOfSprays (map ((skewJackPolTerm n which lambda mu) . fromPartition) nus)
--   where
--     lambda' = toPartitionUnsafe lambda
--     nus = filter (isSuperPartitionOf lambda') (partitions (sum lambda - sum mu))

-- psCombinationsProduct :: 
--   (Eq a, AlgRing.C a) => Map Partition a -> Map Partition a -> Map Partition a 
-- psCombinationsProduct psCombo1 psCombo2 = 
--   DM.fromListWith (AlgAdd.+)
--     [
--       (sortBy (flip compare) (lambda1 ++ lambda2), coeff1 AlgRing.* coeff2)
--       | (lambda1, coeff1) <- DM.assocs psCombo1, 
--         (lambda2, coeff2) <- DM.assocs psCombo2
--     ]