package AGI

import scala.util.Random
import scala.math.{abs,min,floor,sqrt,pow}
import scala.io.Source
import javafx.collections.{ObservableList, FXCollections}
import scala.collection.mutable.{Buffer, ListBuffer}
import scalafx.beans.property.StringProperty

class MDS_Plane (val url1:String,val url2:String,val url3:String){
  def sList2dList(s:List[String],d:List[List[Double]]):List[List[Double]]={
    s match {
      case Nil => List(Nil)
      case (e :: Nil) => d :+ e.split(",").toList.map(x=>x.toDouble)
      case (e :: rst) => sList2dList(rst,d:+e.split(",").toList.map(x=>x.toDouble))
    }
  }
  def sList2iList(s:List[String],d:List[List[Int]]):List[List[Int]]={
    s match {
      case Nil => List(Nil)
      case (e :: Nil) => d :+ e.split(",").toList.map(x=>x.toInt)
      case (e :: rst) => sList2iList(rst,d:+e.split(",").toList.map(x=>x.toInt))
    }
  }

	val p:List[List[Double]] = sList2dList(Source.fromFile(url1).getLines.toList,Nil)
  val adjacency:List[List[Int]] = sList2iList(Source.fromFile(url2).getLines.toList,Nil)
  val eigvals:List[Double] = Source.fromFile(url3).getLines.toList.map(x=>x.toDouble)
  var Xs:ListBuffer[Double] = ListBuffer()
  var Ys:Buffer[Double] = Buffer()
  var boundingX:Double = 0
  var boundingY:Double = 0
		
  val n:Int = p.length
  def countP(lst:List[Double], ans:Int):Int={
    lst match {
      case Nil => ans
      case e::rst => if(e>0) countP(rst,ans+1)
                      else ans
    }
  }
  val d:Int = countP(eigvals,0)
  def takeD(lstlst:List[List[Double]],ans:List[List[Double]],d:Int):List[List[Double]]={
    lstlst match {
      case Nil => ans
      case lst::rst => takeD(rst,ans :+ lst.take(d),d)
    }
  }
  val points:List[List[Double]] = takeD(p,Nil,d)
  val ePos = eigvals.take(d)
  var f1:List[Double] = List()
  var f2:List[Double] = List()
  var flag = true
  for(i <- 0 to d-1){
    if(flag){
      f1 = f1 :+ sqrt(ePos(i))
      f2 = f2 :+ 0d
      flag = false
    } else {
      f2 = f2 :+ sqrt(ePos(i))
      f1 = f1 :+ 0d
      flag = true
    }
  }

  def dot(lst1:List[Double], lst2:List[Double], sum:Double):Double={
    lst1 match {
      case Nil => sum
      case e1::rst1 => lst2 match{
        case Nil => Math.sqrt(sum)
        case e2::rst2 => dot(rst1,rst2, sum + e1*e2)
      }
    }
  }
  var e1 = f1.map(x => x/Math.sqrt(dot(f1,f1,0)))
  var e2 = f2.map(x => x/Math.sqrt(dot(f2,f2,0)))
  var temp1 = e1
  var temp2 = e2

  for(i <- 0 to n-1){
    Xs = Xs :+ dot(points(i),e1,0)
    Ys = Ys :+ dot(points(i),e2,0)
  }

  def boundingBox():Unit={
    this.boundingX = 2 * List(abs(Xs.min),Xs.max).max
    this.boundingY = 2 * List(abs(Ys.min),Ys.max).max
  }

  // 射影平面の更新
  /*def agi(i:int, x2:int, y2:int):(e1,e2)={
    val p = P[i,]
    val x = Xs(i)
    val y = Ys(i)
    val p2 = norm(P[i,])^2
    val m1 = (x*x*y*y)/((p2-y*y)*(p2-y*y)) - (p2 - x*x)
    val m2 = x*y / (p2 - y*y)

    def newton(v0:Double):Double={
      return v0 - f(v0)/f(v0) // 漸化式
    }

    def f(a1:Double):Double={
      return (cul_b1(a1)+y*cul_c1(a1,cul_b1(a1)))*(cul_a2(a1)+x*cul_c2(cul_a2(a1),cul_b2(cul_a2(a1)))-1)
        -(cul_b2(cul_a2(a2))+y*cul_c2(cul_a2(a1),cul_b2(cul_a2(a1))))*(a1 + x*cul_c1(a1,cul_b1(a1)-1))
    }

    def cul_a2(a1:Double):Double={
      return Math.sqrt((y2*y2 - p2)*(x*y*a1 + (p2-y*y)*cul_b1(a1))*(x*y*a1 + (p2-y*y)*cul_b1(a1))
        /((p2 - x*x)*a1 + m2*(x*y*a1+(p2-y*y)*cul_b1(a1)) - x*y*cul_b1(a1))^2 + m1)
    }
    def cul_b1(a1:Double):Double={
      return Math.sqrt(p2 - x2*x2 + m1*a1*a1)+ m2*a1
    }
    def cul_b2(a2:Double):Double={
      return Math.sqrt(p2 - y2*y2 + m1*a2*a2)+ m2*a2
    }
    def cul_c1(a1:Double,b1:Double):Double={
      return (x2 - a1*x - b1*y)/p2
    }
    def cul_c1(a2:Double,b2:Double):Double={
      return (y2 - a2*x - b2*y)/p2
    }
    
    // 制約解消
    val a1 = newton(0)
    val a2 = cul_a2(a1)
    val b1 = cul_b1(a1)
    val b2 = cul_b2(a2)
    val c1 = cul_c1(a1,b1)
    val c2 = cul_c2(a2,b2)

    return (a1*temp1 + b1*temp2 + c1*p, a2*temp1 + b2*temp2 + c2*p)
  }*/

  def update(x2:Double, y2:Double, i:String):Unit={
    val ith = i.toInt
    Xs(ith) = x2
    Ys(ith) = y2
  }

  /*def print():Unit={
    for(i <- 0 to n-1){
      System.out.println(Xs(i))
    }
  }*/
}