package cn.edu.gxust.jiweihuang.scala.math.function.univariate

import cn.edu.gxust.jiweihuang.scala.math.function._
import org.hipparchus.analysis.ParametricUnivariateFunction
import org.hipparchus.analysis.differentiation.DerivativeStructure

import scala.math.{exp, log, pow}

/**
  * <p>reference: https://en.wikipedia.org/wiki/Logistic_function</p>
  * <p>The formula: {{{l(x) = m/(1+exp(-k*(x-x0)))}}}.</p>
  * <p>The property {{{logisticM}}}: The parameter {m} of logistic function.</p>
  * <p>The property {{{logisticK}}}: The parameter {k} of logistic function.</p>
  * <p>The property {{{logisticX0}}}: The parameter {x0} of logistic function.</p>
  *
  * <p>The default {{{logisticM = 1}}},</p>
  * <p>The default {{{logisticK = 1}}}</p>
  * <p>The default {{{logisticX0 = 0}}}</p>
  *
  * <p>require {{{logisticM != 0}}}</p>
  */
class Logistic(val logisticM: Double = 1.0,
               val logisticK: Double = 1.0,
               val logisticX0: Double = 0.0) extends UnivariateFunction
  with UnivariateDifferentiableFunction
  with UnivariateDerivativeFunction
  with UnivariateIntegrableFunction
  with UnivariateIntegralFunction {

  if (logisticM == 0.0) throw new IllegalArgumentException(
    s"Expected the property {logisticM != 0},but get {logisticM = $logisticM}")

  /**
    * The string form of analysis formula of univariate function.
    */
  override val formula: String = s"$logisticM / (1 + exp(-$logisticK * (x - $logisticX0)))"

  def logisticExp(x: Double): Double = {
    Logistic.logisticExp(logisticM, logisticK, logisticX0)(x)
  }

  def logisticExpPluOne(x: Double): Double = {
    Logistic.logisticExpPlusOne(logisticM, logisticK, logisticX0)(x)
  }

  override def derivative(x: Double): Double = {
    Logistic.logisticDerivative(logisticM, logisticK, logisticX0)(x)
  }

  override def integrate(x: Double): Double = {
    Logistic.logisticIntegrate(logisticM, logisticK, logisticX0)(x)
  }

  override def value(t: DerivativeStructure): DerivativeStructure = {
    t.subtract(logisticX0).multiply(-logisticK).exp().add(1).pow(-1).multiply(logisticM)
  }

  override def value(x: Double): Double = {
    Logistic.logistic(logisticM, logisticK, logisticX0)(x)
  }

  override def equals(other: Any): Boolean = other match {
    case that: Logistic =>
      (that canEqual this) &&
        logisticM == that.logisticM &&
        logisticK == that.logisticK &&
        logisticX0 == that.logisticX0
    case _ => false
  }

  override def hashCode(): Int = {
    val state = Seq(logisticM, logisticK, logisticX0)
    state.map(_.hashCode()).foldLeft(0)((a, b) => 31 * a + b)
  }

  def canEqual(other: Any): Boolean = other.isInstanceOf[Logistic]
}

object Logistic {
  /**
    * lexp(x) = exp(-k*(x-x0))
    */
  def logisticExp(logisticM: Double = 1.0,
                  logisticK: Double = 1.0,
                  logisticX0: Double = 0.0)(x: Double): Double = {
    exp(-logisticK * (x - logisticX0))
  }

  /**
    * lexpp1(x) = 1 + exp(-k*(x-x0))
    */
  def logisticExpPlusOne(logisticM: Double = 1.0,
                         logisticK: Double = 1.0,
                         logisticX0: Double = 0.0)(x: Double): Double = {
    1 + logisticExp(logisticM, logisticK, logisticX0)(x)
  }

  /**
    * lexpp1p2(x) = pow(1 + exp(-k*(x-x0)),2)
    */
  def logisticExpPlusOnePow2(logisticM: Double = 1.0,
                             logisticK: Double = 1.0,
                             logisticX0: Double = 0.0)(x: Double): Double = {
    pow(logisticExpPlusOne(logisticM, logisticK, logisticX0)(x), 2)
  }

  /**
    * lexpp1pn1(x) = pow(1 + exp(-k*(x-x0)),-1) = 1/(1 + exp(-k*(x-x0)))
    */
  def logisticExpPlusOnePowN1(logisticM: Double = 1.0,
                              logisticK: Double = 1.0,
                              logisticX0: Double = 0.0)(x: Double): Double = {
    1 / logisticExpPlusOne(logisticM, logisticK, logisticX0)(x)
  }

  /**
    * l(x) = m / (1 + exp(-k*(x-x0)))
    */
  def logistic(logisticM: Double = 1.0,
               logisticK: Double = 1.0,
               logisticX0: Double = 0.0)(x: Double): Double = {
    logisticM / logisticExpPlusOne(logisticM, logisticK, logisticX0)(x)
  }

  /**
    * il(x) = m * (x+log(1 + exp(-k*(x-x0)))/k)
    */
  def logisticIntegrate(logisticM: Double = 1.0,
                        logisticK: Double = 1.0,
                        logisticX0: Double = 0.0)(x: Double): Double = {
    logisticM * (x + log(logisticExpPlusOne(logisticM, logisticK, logisticX0)(x)) / logisticK)
  }

  /**
    * dl(x)
    */
  def logisticDerivative(logisticM: Double = 1.0,
                         logisticK: Double = 1.0,
                         logisticX0: Double = 0.0)(x: Double): Double = {
    (logisticM * logisticK * logisticExp(logisticM, logisticK, logisticX0)(x)) / logisticExpPlusOnePow2(logisticM, logisticK, logisticX0)(x)
  }

  /**
    * dlm(x)
    */
  def logisticDerivativeM(logisticM: Double = 1.0,
                          logisticK: Double = 1.0,
                          logisticX0: Double = 0.0)(x: Double): Double = {
    logisticExpPlusOnePowN1(logisticM, logisticK, logisticX0)(x)
  }

  /**
    * dlk(x)
    */
  def logisticDerivativeK(logisticM: Double = 1.0,
                          logisticK: Double = 1.0,
                          logisticX0: Double = 0.0)(x: Double): Double = {
    -(logisticM * (-x + logisticX0) * logisticExp(logisticM, logisticK, logisticK)(x)) / logisticExpPlusOnePow2(logisticM, logisticK, logisticX0)(x)
  }

  /**
    * dlx0(x)
    */
  def logisticDerivativeX0(logisticM: Double = 1.0,
                           logisticK: Double = 1.0,
                           logisticX0: Double = 0.0)(x: Double): Double = {
    -(logisticK * logisticM * logisticExp(logisticM, logisticK, logisticX0)(x)) / logisticExpPlusOnePow2(logisticM, logisticK, logisticX0)(x)
  }

  class Parametric extends ParametricUnivariateFunction {
    override def value(x: Double, parameters: Double*): Double = {
      checkParameter(parameters: _*)
      logistic(parameters(0), parameters(1), parameters(2))(x)
    }

    override def gradient(x: Double, parameters: Double*): Array[Double] = {
      checkParameter(parameters: _*)
      Array[Double](logisticDerivativeM(parameters(0), parameters(1), parameters(2))(x),
        logisticDerivativeK(parameters(0), parameters(1), parameters(2))(x),
        logisticDerivativeX0(parameters(0), parameters(1), parameters(2))(x))
    }

  }

  def checkParameter(parameters: Double*): Unit = {
    if (parameters == null) throw new IllegalArgumentException(
      s"Expected the parameter {parameters != null},but got {parameters = null}}")
    if (parameters.length != 3) throw new IllegalArgumentException(
      s"Expected the parameter {parameters.length == 3},but got {parameters.length = ${parameters.length}}")
    if (parameters.head == 0) throw new IllegalArgumentException(
      s"Expected the parameter {parameters(0) != 0},but got {parameters(0) = ${parameters.head}}")
  }

  def apply(logisticM: Double = 1.0,
            logisticK: Double = 1.0,
            logisticX0: Double = 0.0): Logistic = {
    Logistic(logisticM, logisticK, logisticX0)
  }

  def unapply(logistic: Logistic): Option[(Double, Double, Double)] =
    if (logistic == null) None
    else Some(logistic.logisticM, logistic.logisticK, logistic.logisticX0)
}
