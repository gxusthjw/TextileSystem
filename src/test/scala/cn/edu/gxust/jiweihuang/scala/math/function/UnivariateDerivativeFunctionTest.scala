package cn.edu.gxust.jiweihuang.scala.math.function

import cn.edu.gxust.jiweihuang.scala.test.UnitSpec

class UnivariateDerivativeFunctionTest extends UnitSpec{
  "UnivariateFunction" should "be right." in {
    import math._
    val udf: UnivariateDerivativeFunction = new UnivariateDerivativeFunction {
      override def formula(): String = "2.0 * pow(3.0 - x, 2.0) + 1.0"
      override def value(x: Double): Double = 2.0 * pow(3.0 - x, 2.0) + 1.0
      override def derivative(x: Double): Double = 4.0 *(3.0 - x)
    }
    println(s"udf.formula() = ${udf.formula()}")
    println(s"udf(3) = ${udf(3)}")
    println(s"udf.value(3) = ${udf.value(3)}")
    println(s"udf.derivative(3) = ${udf.derivative(3)}")
    assertResult(udf(3))(udf.value(3))
  }
}
