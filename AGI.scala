import scala.math.{abs,min}
import scala.collection.mutable.{MutableList, Map}
import scala.util.Random

import breeze.linalg._

import scalafx.Includes._
import scalafx.application.JFXApp
import scalafx.application.JFXApp.PrimaryStage
import scalafx.event.ActionEvent
import scalafx.geometry.Point2D
import scalafx.scene.Scene
import scalafx.scene.control.{Button, ButtonBase, ColorPicker, ToggleButton, ToolBar, ToggleGroup}
import scalafx.scene.effect.DropShadow
import scalafx.scene.input.MouseEvent
import scalafx.scene.layout.{Pane, BorderPane}
import scalafx.scene.paint.Color
import scalafx.scene.shape.{Ellipse, Line, Polygon, Polyline, Rectangle, Shape}

abstract class TDShape
case class TDNoShape()                  extends TDShape
case class Link     (shape:Line)      extends TDShape
case class Node  (shape:Ellipse)   extends TDShape



object AGI extends JFXApp {
  val x, y, x2, y2 = 0
  // constraints
  def c1(a1:Float,b1:Float):Float={
    return a1+b1
  }
  def c2(a2:Float,b2:Float):Float={
    return a2+b2
  }
  def b1(a1:Float):Float={
    return a1+c1(a1,this)
  }
  def b2(a2:Float):Float={}
  def a1(a2:Float):Float={
    return a2+b1(this)
  }

  // Culcurate init (multi-dimensional plane)
  val n:Int = 10  
  val d = DenseMatrix.ones[Double](n,n) - diag(DenseVector.ones[Double](n))

  val A = DenseMatrix.zeros[Double](n,n)
  val dr = DenseVector.zeros[Double](n)
  val dc = DenseVector.zeros[Double](n)
  for(i <- 0 to n-1){
    for(j <- 0 to n-1){  
      dr(i) += d(j,i) * d(j,i)
      dc(i) += d(i,j) * d(i,j) 
    }
  }
  var da2:Double = 0
  for(i <- 0 to n-1){ da2 += dr(i)}
  da2 = da2 / (n*n)
  val dr2 = dr / n.toDouble
  val dc2 = dc / n.toDouble

  for(i <- 0 to n-1){
    for(j <- 0 to n-1){  
      A(i,j) = ((dc2(i)+dr2(j)) - da2 - d(i,j)*d(i,j))/2
    }
  }

  /*
  L,X = np.linalg.eigh(a)
  L.sort()
  L = L[::-1]
  Ln = np.sqrt(np.diag(L))
  val es = eigSym(A)
  val lambda = es._1
  val Ln = diag(lambda)
  val evs = es._2

  P = X * Ln

  f2 = np.array(L)
  f2[::2] = 0
  f1 = L - f2
  e1 = f1 / np.linalg.norm(f1)
  e2 = f2 / np.linalg.norm(f2)*/


  // Update the plane
  /*
  var strokeColor: Color = Color.Black
  var fillColor: Color   = Color.White
  var shapes: MutableList[TDShape] = MutableList()

  type MouseHandler = MouseEvent => Unit

  object mouseAction {
    val _noAction:   MouseHandler = { _ => () }

    var pressHandlers   = Map[String, MouseHandler]()
    var dragHandlers    = Map[String, MouseHandler]()
    var releaseHandlers = Map[String, MouseHandler]()

    def set(id: String) {
      def lookup(handlers: Map[String, MouseHandler]) = {
        try { handlers(id) } catch { case _: NoSuchElementException => _noAction }
      }

      drawingPane.onMousePressed  = lookup(pressHandlers)
      drawingPane.onMouseDragged  = lookup(dragHandlers)
      drawingPane.onMouseReleased = lookup(releaseHandlers)
    }
  }

  object NodeControl { // 楕円の描画に関わるデータ構造と処理（未完成）

    var e = new Ellipse {}
    var p0 = new Point2D(0, 0)

    def onPress(ev: MouseEvent) {
      p0 = new Point2D(ev.x, ev.y)
      e = new Ellipse {
        centerX = p0.x; centerY = p0.y
        stroke = Color.Yellow; fill = Color.Transparent
      }
      drawingPane.content += e
      shapes += Node(e)
    }

    def onDrag(ev: MouseEvent) {
      e.centerX = (p0.x + ev.x) / 2 ; e.centerY = (p0.y + ev.y) / 2
      e.radiusX = abs(p0.x - e.centerX()) ; e.radiusY = abs(p0.y - e.centerY())
    }

    def onRelease(ev: MouseEvent) {
      e.stroke = strokeColor
      e.fill = fillColor
    }
  }

  mouseAction.pressHandlers   += (("Ellipse", NodeControl.onPress))
  mouseAction.dragHandlers    += (("Ellipse", NodeControl.onDrag))
  mouseAction.releaseHandlers += (("Ellipse", NodeControl.onRelease))

  object LinkControl { // 直線の描画に関するデータ構造と処理
  
    var l = new Line {}
    var p0 = new Point2D(0,0)

    def onPress(ev: MouseEvent) {
      p0 = new Point2D(ev.x, ev.y)
      l = new Line {
        startX = p0.x; startY = p0.y
        endX   = p0.x; endY   = p0.y
        stroke = Color.Red; fill = Color.Transparent
        strokeWidth = 3
      }
      drawingPane.content += l
      shapes += Link(l)
    }

    def onDrag(ev: MouseEvent) {
      l.startX = p0.x; l.startY = p0.y
      l.endX   = ev.x; l.endY   = ev.y
    }

    def onRelease(ev: MouseEvent) {
      l.stroke = strokeColor
      l.fill = fillColor
    }

  }

  mouseAction.pressHandlers   += (("Line", LinkControl.onPress))
  mouseAction.dragHandlers    += (("Line", LinkControl.onDrag))
  mouseAction.releaseHandlers += (("Line", LinkControl.onRelease))

  object SelectControl { // 選択ツールに関わるデータ構造と処理
    var selection: TDShape = TDNoShape()
    var p0 = new Point2D(0,0)
    var p1 = new Point2D(0,0)
    var p2 = new Point2D(0,0)

    def onPress(ev: MouseEvent) {
      // 影を消す操作
      selection match {
        case Node(e)   => e.effect = null
        case Link(l)      => l.effect = null
        case TDNoShape()    => ()
      }

      // 影をつける動作
      val x = ev.x; val y = ev.y
      val oShape = (shapes.reverse).find((shape: TDShape) =>
        shape match {
          case Node(e)   => e.contains(x, y)
          case Link(l)      => l.contains(x, y)
        })
      selection = oShape match {
        case Some(shape) => shape
        case _ => TDNoShape()
      }
      selection match {
        case Node(e)   => e.effect = new DropShadow(10, Color.Blue); p1 = new Point2D(e.centerX(), e.centerY())
        case Link(l)      => l.effect = new DropShadow(10, Color.Blue); p1 = new Point2D(l.startX(),l.startY()); p2 = new Point2D(l.endX(),l.endY())
        case TDNoShape()    => ()
      }  
      p0 = new Point2D(ev.x, ev.y)
    }    

    def onDrag(ev: MouseEvent) {
      selection match {
        case Node(e)   => 
          e.centerX = p1.x + ev.x - p0.x ; e.centerY = p1.y + ev.y - p0.y
        case Link(l)      => 
          l.startX = p1.x + ev.x - p0.x ; l.startY = p1.y + ev.y - p0.y
          l.endX   = p2.x + ev.x - p0.x ; l.endY   = p2.y + ev.y - p0.y
        case TDNoShape()    => ()
      }
    }

  }

  mouseAction.pressHandlers   += (("Select", SelectControl.onPress))
  mouseAction.dragHandlers    += (("Select", SelectControl.onDrag))*/

  val drawingPane = new Pane { }

  /*val shapeGroup = new ToggleGroup()

  shapeGroup.selectedToggle.onChange {
    val id = shapeGroup.selectedToggle().asInstanceOf[javafx.scene.control.ToggleButton].id()
    mouseAction.set(id)
  }*/

  /*val shapeTools: List[ToggleButton] = List(
    new ToggleButton {
      id = "Select"; text = "選択"
      graphic = new Rectangle { width = 0; height = 32; fill = Color.Transparent }
      toggleGroup = shapeGroup
      minWidth = 40; minHeight = 40
    },

    new ToggleButton {
      id = "Line"
      graphic = new Line {
        stroke = Color.Black; strokeWidth = 3
        startX = 0; startY = 0; endX = 28;  endY = 28
      }
      toggleGroup = shapeGroup
    }, 
  
    new ToggleButton {
      id = "Ellipse"
      graphic = new Ellipse {
        stroke = Color.Black; fill = Color.White
        radiusX = 16; radiusY = 16
      }
      toggleGroup = shapeGroup
    })*/

  /*val colorTools = Seq(
    new ColorPicker(strokeColor) {
      onAction = { e: ActionEvent => 
        val c = value()
        strokeColor = Color.hsb(c.hue, c.saturation, c.brightness, 0.5)
        SelectControl.selection match {
          case Node(e)   => e.stroke = strokeColor
          case Link(l)      => l.stroke = strokeColor
          case _ => ()
        }
      }
    },
    new ColorPicker(fillColor) {
      onAction = { e: ActionEvent =>
        val c = value()
        fillColor = Color.hsb(c.hue, c.saturation, c.brightness, 0.5)
        SelectControl.selection match {
          case Node(e)   => e.fill = fillColor
          case Link(l)      => l.fill = fillColor
          case _ => ()
        }
      }
    })*/

  stage = new PrimaryStage {
    title = "AGI"
    scene = new Scene(600, 400) {
      root = new BorderPane {
        //top = new ToolBar { content = shapeTools ++ colorTools }
        center = drawingPane
      }
    }
  }
}
