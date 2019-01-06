package cn.edu.gxust.jiweihuang.scala.math

import org.hipparchus.analysis.UnivariateFunction
import org.hipparchus.analysis.differentiation.UnivariateDifferentiableFunction

package object function {

  /**
    * The univariate function
    * <p>
    * it inherit from {{{org.hipparchus.analysis.UnivariateFunction}}}
    * in order to get the functionality of {{{hipparchus}}} library.
    * </p>
    *
    * The follow code is used to usage:
    *
    * {{{
    *   //Create UnivariateFunction instance.
    *   val uf:UnivariateFunction
    *
    *   //calculating the value of function with independent variable x.
    *   val functionValue1:Double = uf(x)
    *
    *   //calculating the value of function with independent variable x in java way.
    *   val functionValue2:Double = uf.value(x)
    *
    *   assert functionValue1 == functionValue2
    * }}}
    **/
  trait TUnivariateFunction extends UnivariateFunction
    with Serializable {

    /**
      * calculating the value of function in object invocation syntax.
      *
      * @param x independent variable.
      * @return the function value
      */
    def apply(x: Double): Double = value(x)

    /**
      * String format of function.
      */
    def formula(): String
  }

  /** The derivative function of univariate function */
  trait TUnivariateDerivativeFunction extends TUnivariateFunction {
    /**
      * calculating the value of derivative function.
      *
      * @param x independent variable.
      * @return the value of derivative function.
      */
    def derivative(x: Double): Double
  }

  /** The integral function of univariate function */
  trait TUnivariateIntegralFunction extends TUnivariateFunction {

    /**
      * calculating the value of integral function.
      * the integral constant is 0.
      *
      * @param x independent variable.
      * @return the value of Integral function.
      */
    def integrate(x: Double): Double

    //definite integration
    /**
      * calculating definite integral value of integral function.
      *
      * @param lowerX lower limit of integral
      * @param upperX upper limit of integral
      * @return definite integral value of integral function.
      */
    def integrate(lowerX: Double, upperX: Double): Double =
      integrate(upperX) - integrate(lowerX)
  }

  /** The inverse function of univariate function */
  trait TUnivariateInverseFunction extends TUnivariateFunction {
    def inverse(x: Double): Array[Double]
  }

  /**
    * The differentiable univariate function
    * The usage:
    * {{{
    * val udf:UnivariateDifferentiableFunction
    * val firstOrderDerivative:Double = udf.differential(x)
    * val secondOrderDerivative:Double = udf.differential(x,2)
    * }}}
    **/
  trait TUnivariateDifferentiableFunction extends TUnivariateFunction
    with UnivariateDifferentiableFunction {

    import org.hipparchus.analysis.differentiation.DSFactory

    def differential(x: Double, order: Int = 1): Double = {
      value(new DSFactory(1, order).variable(0,
        x)).getPartialDerivative(1)
    }
  }

  /** The integrable univariate function */
  trait TUnivariateIntegrableFunction extends TUnivariateFunction {

    import org.hipparchus.analysis.integration._

    private[this] val DefaultIntegrationPointsNumber = 32

    /** Romberg integral algorithm */
    def integrateRomberg(lowerX: Double, upperX: Double, maxIter: Int = Int.MaxValue,
                         relativeAccuracy: Double = BaseAbstractUnivariateIntegrator.DEFAULT_RELATIVE_ACCURACY,
                         absoluteAccuracy: Double = BaseAbstractUnivariateIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                         minimalIterationCount: Int = BaseAbstractUnivariateIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                         maximalIterationCount: Int = RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT
                        ): Double =
      new RombergIntegrator(relativeAccuracy, absoluteAccuracy,
        minimalIterationCount, maximalIterationCount).integrate(maxIter, this, lowerX, upperX)

    /** Simpson integral algorithm */
    def integrateSimpson(lowerX: Double, upperX: Double, maxIter: Int = Int.MaxValue,
                         relativeAccuracy: Double = BaseAbstractUnivariateIntegrator.DEFAULT_RELATIVE_ACCURACY,
                         absoluteAccuracy: Double = BaseAbstractUnivariateIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                         minimalIterationCount: Int = BaseAbstractUnivariateIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                         maximalIterationCount: Int = SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT
                        ): Double =
      new SimpsonIntegrator(relativeAccuracy, absoluteAccuracy,
        minimalIterationCount, maximalIterationCount).integrate(maxIter, this, lowerX, upperX)

    /** MidPoint integral algorithm */
    def integrateMidPoint(lowerX: Double, upperX: Double, maxIter: Int = Int.MaxValue,
                          relativeAccuracy: Double = BaseAbstractUnivariateIntegrator.DEFAULT_RELATIVE_ACCURACY,
                          absoluteAccuracy: Double = BaseAbstractUnivariateIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                          minimalIterationCount: Int = BaseAbstractUnivariateIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                          maximalIterationCount: Int = MidPointIntegrator.MIDPOINT_MAX_ITERATIONS_COUNT
                         ): Double =
      new MidPointIntegrator(relativeAccuracy, absoluteAccuracy,
        minimalIterationCount, maximalIterationCount).integrate(maxIter, this, lowerX, upperX)

    /** Trapezoid integral algorithm */
    def integrateTrapezoid(lowerX: Double, upperX: Double, maxIter: Int = Int.MaxValue,
                           relativeAccuracy: Double = BaseAbstractUnivariateIntegrator.DEFAULT_RELATIVE_ACCURACY,
                           absoluteAccuracy: Double = BaseAbstractUnivariateIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                           minimalIterationCount: Int = BaseAbstractUnivariateIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                           maximalIterationCount: Int = TrapezoidIntegrator.TRAPEZOID_MAX_ITERATIONS_COUNT
                          ): Double =
      new TrapezoidIntegrator(relativeAccuracy, absoluteAccuracy,
        minimalIterationCount, maximalIterationCount).integrate(maxIter, this, lowerX, upperX)

    /** Legendre-Gauss quadrature rule on sub-intervals of the integration range */
    def integrateIterativeLegendreGauss(lowerX: Double, upperX: Double, maxIter: Int = Int.MaxValue,
                                        integrationPointsNumber: Int = DefaultIntegrationPointsNumber,
                                        relativeAccuracy: Double = BaseAbstractUnivariateIntegrator.DEFAULT_RELATIVE_ACCURACY,
                                        absoluteAccuracy: Double = BaseAbstractUnivariateIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                                        minimalIterationCount: Int = BaseAbstractUnivariateIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                                        maximalIterationCount: Int = BaseAbstractUnivariateIntegrator.DEFAULT_MAX_ITERATIONS_COUNT
                                       ): Double =
      new IterativeLegendreGaussIntegrator(integrationPointsNumber,
        relativeAccuracy, absoluteAccuracy, minimalIterationCount,
        maximalIterationCount).integrate(maxIter, this, lowerX, upperX)
  }

}
