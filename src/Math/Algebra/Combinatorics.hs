{-|
Module      : Math.Algebra.Combinatorics
Description : 
Copyright   : (c) Stéphane Laurent, 2024
License     : GPL-3
Maintainer  : laurent_step@outlook.fr

This module provides some functions to compute Kostka numbers with a Jack
parameter, possibly skew, some functions to enumerate semistandard tableaux,
possibly skew, with a given shape and a given weight, and a function to 
enumerate the Gelfand-Tsetlin patterns defined by a skew partition.
-}

module Math.Algebra.Combinatorics
  (
  -- * Kostka numbers
    kostkaNumbers
  , symbolicKostkaNumbers
  , skewKostkaNumbers
  , symbolicSkewKostkaNumbers
  -- * Tableaux
  , semiStandardTableauxWithGivenShapeAndWeight
  , skewTableauxWithGivenShapeAndWeight
  -- * Gelfand-Tsetlin patterns
  , skewGelfandTsetlinPatterns
  ) where
import qualified Data.Foldable                    as DF
import           Data.Map.Strict                  ( 
                                                    Map
                                                  )
import qualified Data.Map.Strict                  as DM
import           Data.Tuple.Extra                 ( 
                                                    second  
                                                  )
import           Math.Algebra.Hspray              (
                                                    RatioOfQSprays
                                                  , unitRatioOfSprays
                                                  )
import           Math.Algebra.Jack.Internal       ( 
                                                    Partition
                                                  , _isPartition
                                                  , _kostkaNumbers
                                                  , _symbolicKostkaNumbers
                                                  , isSkewPartition
                                                  , skewJackInMSPbasis
                                                  , skewSymbolicJackInMSPbasis
                                                  , _skewGelfandTsetlinPatterns
                                                  , _skewTableauxWithGivenShapeAndWeight
                                                  , _semiStandardTableauxWithGivenShapeAndWeight
                                                  )
import           Math.Combinat.Tableaux.Skew      (
                                                    SkewTableau (..)
                                                  )

-- | Kostka numbers \(K_{\lambda,\mu}(\alpha)\) for a given weight of the 
-- partitions \(\lambda\) and \(\mu\) and a given Jack parameter 
-- \(\alpha\) (these are the standard Kostka numbers when
-- \(\alpha=1\)). This returns a map whose keys represent the 
-- partitions \(\lambda\) and the value attached to a partition \(\lambda\)
-- represents the map \(\mu \mapsto K_{\lambda,\mu}(\alpha)\) where the 
-- partition \(\mu\) is included in the keys of this map if and only if 
-- \(K_{\lambda,\mu}(\alpha) \neq 0\).
kostkaNumbers :: 
     Int      -- ^ weight of the partitions
  -> Rational -- ^ Jack parameter
  -> Map Partition (Map Partition Rational)
kostkaNumbers weight alpha 
  | weight < 0 = 
      error "kostkaNumbers: negative weight."
  | weight == 0 =
      DM.singleton [] (DM.singleton [] 1)
  | otherwise =
      _kostkaNumbers weight weight alpha 'P'

-- | Kostka numbers \(K_{\lambda,\mu}(\alpha)\) with symbolic Jack parameter \(\alpha\) 
-- for a given weight of the partitions \(\lambda\) and \(\mu\). This returns a map 
-- whose keys represent the 
-- partitions \(\lambda\) and the value attached to a partition \(\lambda\)
-- represents the map \(\mu \mapsto K_{\lambda,\mu}(\alpha)\) where the 
-- partition \(\mu\) is included in the keys of this map if and only if 
-- \(K_{\lambda,\mu}(\alpha) \neq 0\).
symbolicKostkaNumbers :: 
     Int  -- ^ weight of the partitions
  -> Map Partition (Map Partition RatioOfQSprays)
symbolicKostkaNumbers weight
  | weight < 0 = 
      error "symbolicKostkaNumbers: negative weight."
  | weight == 0 =
      DM.singleton [] (DM.singleton [] unitRatioOfSprays)
  | otherwise =
      _symbolicKostkaNumbers weight weight 'P'

-- | Skew Kostka numbers \(K_{\lambda/\mu, \nu}(\alpha)\) with a given Jack 
-- parameter \(\alpha\) and a given skew partition \(\lambda/\mu\). For \(\alpha=1\)
-- these are the ordinary skew Kostka numbers.
-- This returns a map whose keys represent the partitions \(\nu\).
skewKostkaNumbers ::
     Rational  -- ^ Jack parameter
  -> Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Map Partition Rational
skewKostkaNumbers alpha lambda mu 
  | not (isSkewPartition lambda mu) =
      error "skewKostkaNumbers: invalid skew partition."
  | otherwise = 
      DM.map snd (skewJackInMSPbasis alpha 'P' lambda mu)

-- | Skew Kostka numbers \(K_{\lambda/\mu, \nu}(\alpha)\) with symbolic Jack 
-- parameter \(\alpha\) for a given skew partition \(\lambda/\mu\). 
-- This returns a map whose keys represent the partitions \(\nu\).
symbolicSkewKostkaNumbers ::
     Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Map Partition RatioOfQSprays
symbolicSkewKostkaNumbers lambda mu 
  | not (isSkewPartition lambda mu) =
      error "symbolicSkewKostkaNumbers: invalid skew partition."
  | otherwise = 
      DM.map snd (skewSymbolicJackInMSPbasis 'P' lambda mu)

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
