{-|
Module      : Math.Combinatorics.Kostka
Description : 
Copyright   : (c) Stéphane Laurent, 2024
License     : GPL-3
Maintainer  : laurent_step@outlook.fr

This module provides some functions to compute Kostka-Jack numbers, i.e.
Kostka numbers with a Jack parameter, possibly skew.
-}

module Math.Combinatorics.Kostka
  (
  -- * Kostka numbers
    kostkaNumbersWithGivenLambda
  , kostkaNumbers
  , symbolicKostkaNumbersWithGivenLambda
  , symbolicKostkaNumbers
  , skewKostkaNumbers
  , symbolicSkewKostkaNumbers
  ) where
import           Data.Map.Strict                  ( 
                                                    Map
                                                  )
import qualified Data.Map.Strict                  as DM
import           Math.Algebra.Hspray              (
                                                    RatioOfQSprays
                                                  , unitRatioOfSprays
                                                  )
import           Math.Algebra.Jack.Internal       ( 
                                                    Partition
                                                  , _isPartition
                                                  , _kostkaNumbers
                                                  , _symbolicKostkaNumbers
                                                  , _kostkaNumbersWithGivenLambda
                                                  , _symbolicKostkaNumbersWithGivenLambda
                                                  , isSkewPartition
                                                  , skewJackInMSPbasis
                                                  , skewSymbolicJackInMSPbasis
                                                  )

-- | Kostka numbers \(K_{\lambda,\mu}(\alpha)\) with Jack parameter, or 
-- Kostka-Jack numbers, for a given integer partition \(\lambda\) and a 
-- given Jack parameter \(\alpha\). These are the ordinary Kostka numbers when
-- \(\alpha=1\). The function returns a map whose keys represent the 
-- partitions \(\mu\) and the value attached to a partition \(\mu\) is the
-- Kostka-Jack number \(K_{\lambda,\mu}(\alpha)\). The 
-- partition \(\mu\) is included in the keys of this map if and only if 
-- \(K_{\lambda,\mu}(\alpha) \neq 0\). The Kostka-Jack number 
-- \(K_{\lambda,\mu}(\alpha)\) is the coefficient of the monomial symmetric 
-- polynomial \(m_\mu\) in the expression of the \(P\)-Jack polynomial 
-- \(P_\lambda(\alpha)\) as a linear combination of monomial symmetric 
-- polynomials.
kostkaNumbersWithGivenLambda :: 
     Partition -- ^ the integer partition @lambda@
  -> Rational  -- ^ Jack parameter
  -> Map Partition Rational
kostkaNumbersWithGivenLambda lambda alpha 
  | not (_isPartition lambda) = 
      error "kostkaNumbersWithGivenLambda: invalid integer partition."
  | null lambda =
      DM.singleton [] 1
  | otherwise =
      _kostkaNumbersWithGivenLambda (sum lambda) lambda alpha 'P'

-- | Kostka numbers \(K_{\lambda,\mu}(\alpha)\) with Jack parameter, or 
-- Kostka-Jack numbers, for a given weight of the 
-- partitions \(\lambda\) and \(\mu\) and a given Jack parameter 
-- \(\alpha\). These are the standard Kostka numbers when
-- \(\alpha=1\). The function returns a map whose keys represent the 
-- partitions \(\lambda\) and the value attached to a partition \(\lambda\)
-- represents the map \(\mu \mapsto K_{\lambda,\mu}(\alpha)\) where the 
-- partition \(\mu\) is included in the keys of this map if and only if 
-- \(K_{\lambda,\mu}(\alpha) \neq 0\). The Kostka-Jack number 
-- \(K_{\lambda,\mu}(\alpha)\) is the coefficient of the monomial symmetric 
-- polynomial \(m_\mu\) in the expression of the \(P\)-Jack polynomial 
-- \(P_\lambda(\alpha)\) as a linear combination of monomial symmetric 
-- polynomials.
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

-- | Kostka-Jack numbers \(K_{\lambda,\mu}(\alpha)\) with symbolic Jack 
-- parameter \(\alpha\) for a given weight of the partitions \(\lambda\) 
-- and \(\mu\). This function returns a map whose keys represent the 
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

-- | Kostka-Jack numbers \(K_{\lambda,\mu}(\alpha)\) with symbolic Jack 
-- parameter \(\alpha\) for a given integer partition \(\lambda\). This 
-- function returns a map whose keys represent the 
-- partitions \(\mu\) and the value attached to a partition \(\mu\) is the
-- Kostka-Jack number \(K_{\lambda,\mu}(\alpha)\). The 
-- partition \(\mu\) is included in the keys of this map if and only if 
-- \(K_{\lambda,\mu}(\alpha) \neq 0\). 
symbolicKostkaNumbersWithGivenLambda :: 
     Partition -- ^ the integer partition @lambda@
  -> Map Partition RatioOfQSprays
symbolicKostkaNumbersWithGivenLambda lambda 
  | not (_isPartition lambda) = 
      error "symbolicKostkaNumbersWithGivenLambda: invalid integer partition."
  | null lambda =
      DM.singleton [] unitRatioOfSprays
  | otherwise =
      _symbolicKostkaNumbersWithGivenLambda (sum lambda) lambda 'P'

-- | Skew Kostka numbers \(K_{\lambda/\mu, \nu}(\alpha)\) with a given Jack 
-- parameter \(\alpha\) and a given skew partition \(\lambda/\mu\). For \(\alpha=1\)
-- these are the ordinary skew Kostka numbers.
-- The function returns a map whose keys represent the partitions \(\nu\). 
-- The skew Kostka-Jack number \(K_{\lambda/\mu, \nu}(\alpha)\)
-- is the coefficient of the monomial symmetric 
-- polynomial \(m_\nu\) in the expression of the skew \(P\)-Jack polynomial 
-- \(P_{\lambda/\mu}(\alpha)\) as a linear combination of monomial symmetric 
-- polynomials.
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
-- This function returns a map whose keys represent the partitions \(\nu\).
symbolicSkewKostkaNumbers ::
     Partition -- ^ outer partition of the skew partition
  -> Partition -- ^ inner partition of the skew partition
  -> Map Partition RatioOfQSprays
symbolicSkewKostkaNumbers lambda mu 
  | not (isSkewPartition lambda mu) =
      error "symbolicSkewKostkaNumbers: invalid skew partition."
  | otherwise = 
      DM.map snd (skewSymbolicJackInMSPbasis 'P' lambda mu)

