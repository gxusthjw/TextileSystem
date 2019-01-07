package cn.edu.gxust.jiweihuang.scala

package object test {

  import org.scalatest._

  abstract class UnitSpec extends FlatSpec with Matchers with
    OptionValues with Inside with Inspectors

}
