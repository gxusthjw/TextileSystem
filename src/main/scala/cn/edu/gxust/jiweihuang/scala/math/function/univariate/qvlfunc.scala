package cn.edu.gxust.jiweihuang.scala.math.function.univariate

import cn.edu.gxust.jiweihuang.scala.math.function.{TUnivariateDerivativeFunction, TUnivariateDifferentiableFunction, TUnivariateFunction, TUnivariateIntegrableFunction}
import org.hipparchus.analysis.ParametricUnivariateFunction
import org.hipparchus.analysis.differentiation.DerivativeStructure

import scala.math._

trait TQuadraticVertexLogistic extends TUnivariateFunction
  with TUnivariateDifferentiableFunction
  with TUnivariateIntegrableFunction
  with TUnivariateDerivativeFunction {

  val quadraticVertexA: Double = 1
  val quadraticVertexB: Double = 0
  val quadraticVertexC: Double = 0
  val logisticM: Double = 1
  val logisticK: Double = 1
  val logisticX0: Double = 0
  val quadraticVertexLogisticD: Double = 0

  override val formula: String = s"($quadraticVertexA * pow(x - $quadraticVertexB,2) + $quadraticVertexC) * $logisticM / (1 + exp(-$logisticK*(x - $logisticX0))) + $quadraticVertexLogisticD"
  val quadraticVertex: TQuadraticVertex = QuadraticVertex(quadraticVertexA, quadraticVertexB, quadraticVertexC)
  val logistic: TLogistic = Logistic(logisticM, logisticK, logisticX0)

  override def derivative(x: Double): Double = quadraticVertex.derivative(x) * logistic.value(x) + logistic.derivative(x) * quadraticVertex.value(x)

  override def value(ds: DerivativeStructure): DerivativeStructure = quadraticVertex.value(ds).multiply(logistic.value(ds)).add(quadraticVertexLogisticD)

  override def value(x: Double): Double = logistic.value(x) * quadraticVertex.value(x) + quadraticVertexLogisticD
}

object TQuadraticVertexLogistic {

  import cn.edu.gxust.jiweihuang.scala.math.function.univariate.TLogistic._
  import cn.edu.gxust.jiweihuang.scala.math.function.univariate.TQuadraticVertex._

  def quadraticVertexLogistic(a: Double, b: Double, c: Double,
                              m: Double, k: Double, x0: Double,
                              d: Double)(x: Double): Double = {
    quadraticVertex(a, b, c)(x) * logistic(m, k, x0)(x) + d
  }

  def quadraticVertexLogisticDerivative(a: Double, b: Double, c: Double,
                                        m: Double, k: Double, x0: Double,
                                        d: Double)(x: Double): Double = {
    quadraticVertexDerivative(a, b, c)(x) * logistic(m, k, x0)(x) + logisticDerivative(m, k, x0)(x) * quadraticVertex(a, b, c)(x)
  }

  def quadraticVertexLogisticDerivativeA(a: Double, b: Double, c: Double,
                                         m: Double, k: Double, x0: Double,
                                         d: Double)(x: Double): Double = {
    pow(x - b, 2) * logistic(m, k, x0)(x)
  }

  def quadraticVertexLogisticDerivativeB(a: Double, b: Double, c: Double,
                                         m: Double, k: Double, x0: Double,
                                         d: Double)(x: Double): Double = {
    quadraticVertexDerivative(a, b, c)(x) * logistic(m, k, x0)(x)
  }

  def quadraticVertexLogisticDerivativeC(a: Double, b: Double, c: Double,
                                         m: Double, k: Double, x0: Double,
                                         d: Double)(x: Double): Double = {
    logistic(m, k, x0)(x)
  }

  def quadraticVertexLogisticDerivativeM(a: Double, b: Double, c: Double,
                                         m: Double, k: Double, x0: Double,
                                         d: Double)(x: Double): Double = {
    quadraticVertex(a, b, c)(x) / logisticExpPlusOne(m, k, x0)(x)
  }

  def quadraticVertexLogisticDerivativeK(a: Double, b: Double, c: Double,
                                         m: Double, k: Double, x0: Double,
                                         d: Double)(x: Double): Double = {
    logisticExp(k, x0)(x) * quadraticVertex(a, b, c)(x) * (x - x0) * m / logisticExpPlusOnePow2(m, k, x0)(x)
  }

  def quadraticVertexLogisticDerivativeX0(a: Double, b: Double, c: Double,
                                          m: Double, k: Double, x0: Double,
                                          d: Double)(x: Double): Double = {
    -logisticExp(k, x0)(x) * quadraticVertex(a, b, c)(x) * k * m / logisticExpPlusOnePow2(m, k, x0)(x)
  }

  def quadraticVertexLogisticDerivativeD(a: Double, b: Double, c: Double,
                                         m: Double, k: Double, x0: Double,
                                         d: Double)(x: Double): Double = 1

  final class Parametric extends ParametricUnivariateFunction {
    override def value(x: Double, parameters: Double*): Double = {
      checkParameter(parameters: _*)
      TQuadraticVertexLogistic.quadraticVertexLogistic(parameters(0),
        parameters(1), parameters(2), parameters(3), parameters(4),
        parameters(5), parameters(6))(x)
    }

    override def gradient(x: Double, parameters: Double*): Array[Double] = {
      checkParameter(parameters: _*)
      Array[Double](quadraticVertexLogisticDerivativeA(parameters(0),
        parameters(1), parameters(2), parameters(3), parameters(4),
        parameters(5), parameters(6))(x),
        quadraticVertexLogisticDerivativeB(parameters(0),
          parameters(1), parameters(2), parameters(3), parameters(4),
          parameters(5), parameters(6))(x),
        quadraticVertexLogisticDerivativeC(parameters(0),
          parameters(1), parameters(2), parameters(3), parameters(4),
          parameters(5), parameters(6))(x),
        quadraticVertexLogisticDerivativeM(parameters(0),
          parameters(1), parameters(2), parameters(3), parameters(4),
          parameters(5), parameters(6))(x),
        quadraticVertexLogisticDerivativeK(parameters(0),
          parameters(1), parameters(2), parameters(3), parameters(4),
          parameters(5), parameters(6))(x),
        quadraticVertexLogisticDerivativeX0(parameters(0),
          parameters(1), parameters(2), parameters(3), parameters(4),
          parameters(5), parameters(6))(x),
        quadraticVertexLogisticDerivativeD(parameters(0),
          parameters(1), parameters(2), parameters(3), parameters(4),
          parameters(5), parameters(6))(x))
    }

    def checkParameter(parameters: Double*): Unit = {
      if (parameters == null) throw new IllegalArgumentException(
        s"Expected the parameter {parameters != null},but got {parameters = null}}")
      if (parameters.length != 7) throw new IllegalArgumentException(
        s"Expected the parameter {parameters.length == 7},but got {parameters.length = ${parameters.length}}")
      if (parameters.head == 0) throw new IllegalArgumentException(
        s"Expected the parameter {parameters(0) != 0},but got {parameters(0) = ${parameters.head}}")
      if (parameters(3) == 0) throw new IllegalArgumentException(
        s"Expected the parameter {parameters(3) != 0},but got {parameters(3) = ${parameters(3)}}")
    }
  }

  def apply(quadraticVertexA: Double = 1,
            quadraticVertexB: Double = 0,
            quadraticVertexC: Double = 0,
            logisticM: Double = 1,
            logisticK: Double = 1,
            logisticX0: Double = 0,
            quadraticVertexLogisticD: Double = 0): TQuadraticVertexLogistic =
    QuadraticVertexLogistic(quadraticVertexA, quadraticVertexB, quadraticVertexC, logisticM,
      logisticK, logisticX0, quadraticVertexLogisticD)

  def unapply(quadraticVertexLogistic: TQuadraticVertexLogistic): Option[(Double, Double, Double, Double, Double, Double, Double)] =
    if (quadraticVertexLogistic == null) None
    else Some(quadraticVertexLogistic.quadraticVertexA,
      quadraticVertexLogistic.quadraticVertexB,
      quadraticVertexLogistic.quadraticVertexC,
      quadraticVertexLogistic.logisticM,
      quadraticVertexLogistic.logisticK,
      quadraticVertexLogistic.logisticX0,
      quadraticVertexLogistic.quadraticVertexLogisticD)
}

final case class QuadraticVertexLogistic(override val quadraticVertexA: Double = 1.0,
                                         override val quadraticVertexB: Double = 0.0,
                                         override val quadraticVertexC: Double = 0.0,
                                         override val logisticM: Double = 1.0,
                                         override val logisticK: Double = 1.0,
                                         override val logisticX0: Double = 0.0,
                                         override val quadraticVertexLogisticD: Double = 0.0)
  extends TQuadraticVertexLogistic