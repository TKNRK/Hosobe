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

	val P:List[List[Double]] = sList2dList(Source.fromFile(url1).getLines.toList,Nil)
  val adjacency:List[List[Int]] = sList2iList(Source.fromFile(url2).getLines.toList,Nil)
  val eigvals:List[Double] = Source.fromFile(url3).getLines.toList.map(x=>x.toDouble)
  var Xs:ListBuffer[Double] = ListBuffer()
  var Ys:Buffer[Double] = Buffer()
  var boundingX:Double = 0
  var boundingY:Double = 0
		
  val n:Int = P.length
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
  val points:List[List[Double]] = takeD(P,Nil,d)
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

  for(i <- 0 to n-1){
    Xs = Xs :+ dot(points(i),e1,0)
    Ys = Ys :+ dot(points(i),e2,0)
  }

  def boundingBox():Unit={
    this.boundingX = 2 * List(abs(Xs.min),Xs.max).max
    this.boundingY = 2 * List(abs(Ys.min),Ys.max).max
  }

  def update(x2:Double, y2:Double, i:String):Unit={
    val ith = i.toInt
    Xs(ith) = x2
    Ys(ith) = y2
    //val e = this.agi(ith,x2,y2)
    //e1 = e._1 ; e2 = e._2
    boundingBox()
  }

  /*def print():Unit={
    for(i <- 0 to n-1){
      System.out.println(Xs(i))
    }
  }*/

  // 射影平面の更新(newton法)
  def agi(i:int, x2:int, y2:int):(e1,e2)={
    val p = P(i)
    val x = Xs(i)
    val y = Ys(i)
    val p2 = dot(p,p)
    val m1 = (x*x*y*y)/((p2-y*y)*(p2-y*y)) - (p2 - x*x)
    val m2 = x*y / (p2 - y*y)

    var a1:Double = 0d
    var a2:Double = 0d
    var temp_a1:Double = a1
    var temp_a2:Double = a2
    val epsilon:Double = 0.001

    def newton2():Unit={
      temp_a1 = a1 ; temp_a2 = a2
      a1 = temp_a1 - (dg_da2(temp_a1,temp_a2)*f(temp_a1,temp_a2) - df_da2(temp_a1,temp_a2)*g(temp_a1,temp_a2))/delta(temp_a1,temp_a2) // a1 update
      a2 = temp_a2 - (df_da1(temp_a1,temp_a2)*g(temp_a1,temp_a2) - dg_da1(temp_a1,temp_a2)*f(temp_a1,temp_a2))/delta(temp_a1,temp_a2) // a2 update
    }
    def delta(a1:Double,a2:Double):Double={
      return df_da1(a1,a2)*dg_da2(a1,a2) - df_da2(a1,a2)*dg_da1(a1,a2)
    }

    def f(a1:Double,a2:Double):Double={
      val b1 = cul_b1(a1) ; val b2 = cul_b2(a2)
      val c1 = cul_c1(a1,b2) ; val c2 = cul_c2(a2,b2)
      return a1*a2 + b1*b2 + c1*c2*p2 + (a1*c2 + c1*a2)*x + (b1*c2 + c1*b2)*y
    }
    def df_da1(a1:Double,a2:Double):Double={
      val b1 = cul_b1(a1) ; val b2 = cul_b2(a2)
      val c1 = cul_c1(a1,b2) ; val c2 = cul_c2(a2,b2)
      val db1 = cul_db1(a1) ; val dc1 = cul_dc1(db1)
      return a2 + db1*b2 + dc1*c2*p2 + (c2 + dc1*a2)*x + (db1*c2 + dc1*b2)*y
    }
    def df_da2(a1:Double,a2:Double):Double={
      val b1 = cul_b1(a1) ; val b2 = cul_b2(a2)
      val c1 = cul_c1(a1,b2) ; val c2 = cul_c2(a2,b2)
      val db2 = cul_db2(a2) ; val dc2 = cul_dc2(db2)
      return a1 + b1*db2 + c1*dc2*p2 + (a1*dc2 + c1)*x + (b1*dc2 + c1*db2)*y
    }
    def g(a1:Double,a2:Double):Double={
      val b1 = cul_b1(a1) ; val b2 = cul_b2(a2)
      val c1 = cul_c1(a1,b2) ; val c2 = cul_c2(a2,b2)
      return (b1 + c1*y)*(a2 + c2*x - 1) - (a1 + c1*x - 1)*(b2 + c2*y)
    }
    def dg_da1(a1:Double,a2:Double):Double={
      val b1 = cul_b1(a1) ; val b2 = cul_b2(a2)
      val c1 = cul_c1(a1,b2) ; val c2 = cul_c2(a2,b2)
      val db1 = cul_db1(a1) ; val dc1 = cul_dc1(db1)
      return (db1 + dc1*y)*(a2 + c2*x - 1) - dc1*x*(b2 + c2*y)
    }
    def dg_da2(a1:Double,a2:Double):Double={
      val b1 = cul_b1(a1) ; val b2 = cul_b2(a2)
      val c1 = cul_c1(a1,b2) ; val c2 = cul_c2(a2,b2)
      val db2 = cul_db2(a2) ; val dc2 = cul_dc2(db2)
      return (b1 + c1*y)*dc2*x - (a1 + c1*x - 1)*(db2 + dc2*y)
    }

    def cul_b1(a1:Double):Double={
      return Math.sqrt(p2 - x2*x2 + m1*a1*a1)+ m2*a1
    }
    def cul_db1(a1:Double):Double={
      return (-1*(2*m1*a1))/Math.sqrt(p2 - x2*x2 + m1*a1*a1) + m2
    }
    def cul_b2(a2:Double):Double={
      return Math.sqrt(p2 - y2*y2 + m1*a2*a2)+ m2*a2
    }
    def cul_db2(a2:Double):Double={
      return (-1*(2*m1*a2))/Math.sqrt(p2 - y2*y2 + m1*a2*a2) + m2
    }
    def cul_c1(a1:Double,b1:Double):Double={
      return (x2 - a1*x - b1*y)/p2
    }
    def cul_dc1(db1:Double):Double={
      return (-x - db1*y)/p2
    }
    def cul_c2(a2:Double,b2:Double):Double={
      return (y2 - a2*x - b2*y)/p2
    }
    def cul_dc2(db2:Double):Double={
      return -1*(x + db2*y)/p2
    }

    def branch():Boolean={
      (abs(a1 - temp_a1) + abs(a2 - temp_a2)) < epsilon
    }

    // 制約解消
    newton2()    
    val b1 = cul_b1(a1)
    val b2 = cul_b2(a2)
    val c1 = cul_c1(a1,b1)
    val c2 = cul_c2(a2,b2)

    def vec(a:Double,b:Double,c:Double):List[Double]={
      var v:List[Double] = Nil
      for(i <- 0 to n-1){}
      return v
    }

    //return (vec(a1,b1,c1),vec(a1,b1,c1))
    return (p,p)
  }
}