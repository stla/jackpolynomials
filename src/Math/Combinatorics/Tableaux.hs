{-|
Module      : Math.Combinatorics.Tableaux
Description : 
Copyright   : (c) Stéphane Laurent, 2024
License     : GPL-3
Maintainer  : laurent_step@outlook.fr

This module provides some functions to enumerate semistandard tableaux,
possibly skew, with a given shape and a given weight, and a function to 
enumerate the Gelfand-Tsetlin patterns defined by a skew partition.
-}

module Math.Combinatorics.Tableaux
  (
  -- * Tableaux
    semiStandardTableauxWithGivenShapeAndWeight
  , skewTableauxWithGivenShapeAndWeight
  -- * Gelfand-Tsetlin patterns
  , skewGelfandTsetlinPatterns
  ) where
import qualified Data.Foldable                    as DF
import           Data.Tuple.Extra                 ( 
                                                    second  
                                                  )
import           Math.Algebra.Jack.Internal       ( 
                                                    Partition
                                                  , _isPartition
                                                  , isSkewPartition
                                                  , _skewGelfandTsetlinPatterns
                                                  , _skewTableauxWithGivenShapeAndWeight
                                                  , _semiStandardTableauxWithGivenShapeAndWeight
                                                  )
import           Math.Combinat.Tableaux.Skew      (
                                                    SkewTableau (..)
                                                  )

-- | Skew Gelfand-Tsetlin patterns defined by a skew partition and a weight vector.
skewGelfandTsetlinPatterns :: 
     Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> [Int]     -- ^ weight
  -> [[Partition]]
skewGelfandTsetlinPatterns lambda mu weight 
  | not (isSkewPartition lambda mu) =
     error "skewGelfandTsetlinPatterns: invalid skew partition."
  | otherwise = 
      map (map DF.toList) (_skewGelfandTsetlinPatterns lambda mu weight)

-- | Skew semistandard tableaux with a given shape (a skew partition) and
-- a given weight vector. The weight is the vector whose @i@-th element is the 
-- number of occurrences of @i@ in the tableau.
skewTableauxWithGivenShapeAndWeight :: 
     Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> [Int]     -- ^ weight
  -> [SkewTableau Int] 
skewTableauxWithGivenShapeAndWeight lambda mu weight 
  | not (isSkewPartition lambda mu) =
     error "skewTableauxWithGivenShapeAndWeight: invalid skew partition."
  | otherwise = 
      map (SkewTableau . (map (second DF.toList)))
        (_skewTableauxWithGivenShapeAndWeight lambda mu weight)

-- | Semistandard tableaux with a given shape (an integer partition) and
-- a given weight vector. The weight is the vector whose @i@-th element is the 
-- number of occurrences of @i@ in the tableau.
semiStandardTableauxWithGivenShapeAndWeight :: 
     Partition   -- ^ shape, integer partition
  -> [Int]       -- ^ weight
  -> [[[Int]]]
semiStandardTableauxWithGivenShapeAndWeight lambda weight 
  | not (_isPartition lambda) =
      error "semiStandardTableauxWithGivenShapeAndWeight: invalid partition."
  | any (< 0) weight =
      []
  | otherwise = 
      map (map DF.toList) 
        (_semiStandardTableauxWithGivenShapeAndWeight lambda weight)
