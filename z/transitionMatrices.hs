import qualified Data.Map.Strict                as DM
import Math.Combinat.Partitions.Integer         ( 
                                                  toPartition
                                                , fromPartition
                                                , mkPartition
                                                , partitions 
                                                , dualPartition
                                                )
import qualified Math.Combinat.Partitions.Integer as PI
import Math.Combinat.Tableaux.GelfandTsetlin    ( kostkaNumber )
import qualified Math.Combinat.Tableaux.GelfandTsetlin as GT


b_lambda_mu :: [Int] -> [Int] -> Int
b_lambda_mu lambda mu = sum $ DM.elems wholeMap 
  where
    parts = partitions (sum lambda)
    zeros = DM.fromList (zip parts (repeat 0))
    map1 = DM.union (GT.kostkaNumbersWithGivenMu (mkPartition lambda)) zeros
    map2 = DM.union (GT.kostkaNumbersWithGivenMu (mkPartition mu)) zeros
    wholeMap = DM.unionWithKey (\part kn1 _ -> kn1 * (map2 DM.! (dualPartition part))) map1 map2

b_lambda_mu' :: [Int] -> [Int] -> Int
b_lambda_mu' lambda mu = sum $ DM.elems wholeMap 
  where
    parts = partitions (sum lambda)
    map1 = GT.kostkaNumbersWithGivenMu (mkPartition lambda)
    map2 = DM.mapKeys dualPartition (GT.kostkaNumbersWithGivenMu (mkPartition mu)) 
    wholeMap = DM.intersectionWith (*) map1 map2
