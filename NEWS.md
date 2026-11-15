# pacfwei 0.3.0

- Stabilize multivariate calculations using SVD whitening of the centered
  observations. Both implementations transform the residual covariance and
  cross-covariance matrices back before normalization, preserving the original
  componentwise partial autocorrelations.
- Replace general LU solves with Cholesky and triangular solves. Add numerical
  rank, conditioning, and joint residual-covariance validity checks, including
  the final requested lag. Preserve the direct lag-one calculation for
  redundant series with positive marginal variances.
- Use power-of-two prescaling before centering to avoid loss of information
  for observations with large offsets, while retaining numerical-range
  protection for huge and very small finite values.
- Clamp only roundoff-sized correlation-bound excesses; report larger
  violations as errors.
- Apply the same lag validation to univariate and multivariate inputs. Accept
  `lag.max = Inf`, reject malformed arguments consistently, and give the default
  maximum lag a minimum of one for short, wide data.
- Reject abbreviated names such as `meth` for `method`, and check numeric input
  type before observation count.
- Add readable plain-text alternatives for mathematical help expressions and
  document the stabilization, numerical policies, and API changes.
