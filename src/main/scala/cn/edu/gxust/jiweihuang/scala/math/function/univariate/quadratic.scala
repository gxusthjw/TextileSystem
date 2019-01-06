package cn.edu.gxust.jiweihuang.scala.math.function.univariate

import cn.edu.gxust.jiweihuang.scala.math.function._
import org.hipparchus.analysis.ParametricUnivariateFunction
import org.hipparchus.analysis.differentiation.DerivativeStructure

import scala.math._

/**
  * The trait [[TQuadratic]] is used for
  * representing quadratic function.
  *
  * @see TUnivariateFunction
  * @see TUnivariateInverseFunction
  * @see TUnivariateDerivativeFunction
  * @see TUnivariateDifferentiableFunction
  * @see TUnivariateIntegrableFunction
  * @see TUnivariateIntegrateFunction
  */
trait TQuadratic extends TUnivariateFunction
  with TUnivariateInverseFunction
  with TUnivariateDerivativeFunction
  with TUnivariateDifferentiableFunction
  with TUnivariateIntegralFunction
  with TUnivariateIntegrableFunction {

  /**
    * The vertex coordinates (x,y) of quadratic function
    */
  val vertex: Array[Double]

  /**
    * <p>whether quadratic function is invert.</p>
    * <p>if invert (i.e. open side down) return true else return false.</p>
    */
  val isInvert: Boolean

  /**
    * the x coordinate of intersection with x axis of function
    */
  val xIntersection: Array[Double]

  /**
    * the y coordinate of intersection with y axis of function
    */
  val yIntersection: Double

  /**
    * The number of intersection with x axis of function
    */
  val xIntersectionNum: Int
}

/**
  * <p>reference: https://en.wikipedia.org/wiki/Quadratic_function</p>
  * <p>The formula:{{{q(x)=a*(x-b)^2+c}}}.</p>
  *
  * <p>The property: {{{quadraticVertexA}}}:the parameter {a} of the vertex form of quadratic function.</p>
  * <p>The property: {{{quadraticVertexB}}}:the parameter {b} of the vertex form of quadratic function.</p>
  * <p>The property: {{{quadraticVertexC}}}:the parameter {c} of the vertex form of quadratic function.</p>
  *
  * <p>require {{{quadraticVertexA != 0}}}</p>
  *
  * @see TQuadratic
  */
trait TQuadraticVertex extends TQuadratic {
  val quadraticVertexA: Double = 1.0
  val quadraticVertexB: Double = 0.0
  val quadraticVertexC: Double = 0.0
  /**
    * Ensure the parameters {{{quadraticVertexA != 0}}}
    */
  if (quadraticVertexA == 0) throw new IllegalArgumentException(
    s"Expected the parameter {quadraticVertexA != 0},but got {quadraticVertexA = $quadraticVertexA}.")

  override val vertex: Array[Double] = Array(quadraticVertexB, quadraticVertexC)
  /**
    * <p>whether quadratic function is invert.</p>
    * <p>if invert (i.e. open side down) return true else return false.</p>
    */
  override val isInvert: Boolean = quadraticVertexA < 0
  /**
    * the x coordinate of intersection with x axis of function
    */
  override val xIntersection: Array[Double] = {
    val tem = -quadraticVertexC / quadraticVertexA
    if (tem > 0) Array[Double](quadraticVertexB - sqrt(tem), quadraticVertexB + sqrt(tem))
    else if (tem == 0) Array[Double](quadraticVertexB) else Array[Double]()
  }
  /**
    * the y coordinate of intersection with y axis of function
    */
  override val yIntersection: Double = quadraticVertexA * pow(quadraticVertexB, 2) +
    quadraticVertexC

  /**
    * The number of intersection with x axis of function
    */
  override val xIntersectionNum: Int = xIntersection.length
  /**
    * The String of vertex of quadratic function.
    */
  override val formula: String = s"$quadraticVertexA * pow(x - $quadraticVertexB, 2) + $quadraticVertexC"

  override def value(x: Double): Double = {
    TQuadraticVertex.quadraticVertex(quadraticVertexA, quadraticVertexB, quadraticVertexC)(x)
  }

  override def value(x: DerivativeStructure): DerivativeStructure = {
    x.subtract(quadraticVertexB).pow(2).multiply(quadraticVertexA).add(quadraticVertexC)
  }

  override def inverse(y: Double): Array[Double] = {
    TQuadraticVertex.quadraticVertexInverse(quadraticVertexA, quadraticVertexB, quadraticVertexC)(y)
  }

  override def integrate(x: Double): Double = {
    TQuadraticVertex.quadraticVertexIntegrate(quadraticVertexA, quadraticVertexB, quadraticVertexC)(x)
  }

  override def derivative(x: Double): Double = {
    TQuadraticVertex.quadraticVertexDerivative(quadraticVertexA, quadraticVertexB, quadraticVertexC)(x)
  }

  override def equals(other: Any): Boolean = other match {
    case that: TQuadraticVertex =>
      (that canEqual this) &&
        quadraticVertexA == that.quadraticVertexA &&
        quadraticVertexB == that.quadraticVertexB &&
        quadraticVertexC == that.quadraticVertexC
    case _ => false
  }

  def canEqual(other: Any): Boolean = other.isInstanceOf[TQuadraticVertex]

  override def hashCode(): Int = {
    val state = Seq(quadraticVertexA, quadraticVertexB, quadraticVertexC)
    state.map(_.hashCode()).foldLeft(0)((a, b) => 31 * a + b)
  }
}


object TQuadraticVertex {

  /**
    * q(x) = a * pow(x-b,2) + c
    */
  def quadraticVertex(quadraticVertexA: Double = 1.0,
                      quadraticVertexB: Double = 0.0,
                      quadraticVertexC: Double = 0.0)(x: Double): Double = {
    quadraticVertexA * pow(x - quadraticVertexB, 2) + quadraticVertexC
  }

  /**
    * qi1(x) = b - sqrt(a*(y-c))/a
    * qi2(x) = b + sqrt(a*(y-c))/a
    */
  def quadraticVertexInverse(quadraticVertexA: Double = 1.0,
                             quadraticVertexB: Double = 0.0,
                             quadraticVertexC: Double = 0.0)(y: Double): Array[Double] = {
    val tem = quadraticVertexA * (y - quadraticVertexC)
    if (tem < 0) {
      new Array[Double](0)
    } else if (tem == 0) {
      Array[Double](quadraticVertexB)
    } else {
      Array[Double](quadraticVertexB - sqrt(tem) / quadraticVertexA, quadraticVertexB + sqrt(tem) / quadraticVertexA)
    }
  }

  /**
    * dq(x) = 2 * a * (x-b)
    */
  def quadraticVertexDerivative(quadraticVertexA: Double = 1.0,
                                quadraticVertexB: Double = 0.0,
                                quadraticVertexC: Double = 0.0)(x: Double): Double = {
    2 * quadraticVertexA * (x - quadraticVertexB)
  }

  /**
    * iq(x) = a * pow(b,2) * x + c * x - a * b * pow(x,2)  +  (a * pow(x,3)) / 3
    */
  def quadraticVertexIntegrate(quadraticVertexA: Double = 1.0,
                               quadraticVertexB: Double = 0.0,
                               quadraticVertexC: Double = 0.0)(x: Double): Double = {
    quadraticVertexA * pow(quadraticVertexB, 2) * x + quadraticVertexC * x -
      quadraticVertexA * quadraticVertexB * pow(x, 2) +
      quadraticVertexA * pow(x, 3) / 3
  }

  /**
    * dqa(x) = pow(x-b,2)
    */
  def quadraticVertexDerivativeA(quadraticVertexA: Double = 1.0,
                                 quadraticVertexB: Double = 0.0,
                                 quadraticVertexC: Double = 0.0)(x: Double): Double = {
    pow(x - quadraticVertexB, 2)
  }

  /**
    * dqb(x) = -2 * a * (x-b)
    */
  def quadraticVertexDerivativeB(quadraticVertexA: Double = 1.0,
                                 quadraticVertexB: Double = 0.0,
                                 quadraticVertexC: Double = 0.0)(x: Double): Double = {
    -2 * quadraticVertexA * (x - quadraticVertexB)
  }

  /**
    * dqc(x) = 1
    */
  def quadraticVertexDerivativeC(quadraticVertexA: Double = 1.0,
                                 quadraticVertexB: Double = 0.0,
                                 quadraticVertexC: Double = 0.0)(x: Double): Double = 1.0

  def checkParameter(parameters: Double*): Unit = {
    if (parameters == null) throw new IllegalArgumentException(
      s"Expected the parameter {parameters != null},but got {parameters = null}}")
    if (parameters.length != 3) throw new IllegalArgumentException(
      s"Expected the parameter {parameters.length == 3},but got {parameters.length = ${parameters.length}}")
    if (parameters.head == 0) throw new IllegalArgumentException(
      s"Expected the parameter {parameters(0) != 0},but got {parameters(0) = ${parameters.head}}")
  }

  final class Parametric extends ParametricUnivariateFunction {
    override def value(x: Double, parameters: Double*): Double = {
      checkParameter(parameters: _*)
      quadraticVertex(parameters(0), parameters(1), parameters(2))(x)
    }

    override def gradient(x: Double, parameters: Double*): Array[Double] = {
      checkParameter(parameters: _*)
      Array[Double](quadraticVertexDerivativeA(parameters(0), parameters(1), parameters(2))(x),
        quadraticVertexDerivativeB(parameters(0), parameters(1), parameters(2))(x),
        quadraticVertexDerivativeC(parameters(0), parameters(1), parameters(2))(x))
    }

  }

  def apply(quadraticVertexA: Double = 1.0, quadraticVertexB: Double = 0.0,
            quadraticVertexC: Double = 0.0): TQuadraticVertex =
    QuadraticVertex(quadraticVertexA, quadraticVertexB, quadraticVertexC)

  def unapply(qv: TQuadraticVertex): Option[(Double, Double, Double)] = {
    if (qv == null) {
      None
    } else {
      Some(qv.quadraticVertexA, qv.quadraticVertexB, qv.quadraticVertexC)
    }
  }

}

final case class QuadraticVertex(override val quadraticVertexA: Double = 1.0,
                                 override val quadraticVertexB: Double = 0.0,
                                 override val quadraticVertexC: Double = 0.0)
  extends TQuadraticVertex