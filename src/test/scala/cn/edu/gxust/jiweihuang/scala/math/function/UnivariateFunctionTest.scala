package cn.edu.gxust.jiweihuang.scala.math.function

import cn.edu.gxust.jiweihuang.scala.test.UnitSpec

class UnivariateFunctionTest extends UnitSpec {
  "UnivariateFunction" should "be right." in {
    import math._
    val uf: UnivariateFunction = new UnivariateFunction {
      override def formula(): String = "2.0 * pow(3.0 - x, 2.0) + 1.0"
      override def value(x: Double): Double = 2.0 * pow(3.0 - x, 2.0) + 1.0
    }
    println(s"uf.formula() = ${uf.formula()}")
    println(s"uf(3) = ${uf(3)}")
    println(s"uf.value(3) = ${uf.value(3)}")
    assertResult(uf(3))(uf.value(3))
  }
}
