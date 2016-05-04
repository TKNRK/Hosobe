package AGI

import scala.math.{abs,min,floor,sqrt,pow}
import scala.collection.mutable.{Buffer,MutableList, Map}
import scala.compat.Platform
import javafx.collections.{ObservableList, FXCollections}

import scalafx.Includes._
import scalafx.application.JFXApp
import scalafx.application.JFXApp.PrimaryStage
import scalafx.event.ActionEvent
import scalafx.geometry.Point2D
import scalafx.scene.Scene
import scalafx.scene.control.{ListView, Button, ButtonBase, ToggleButton, ToolBar, ToggleGroup}
import scalafx.scene.effect.DropShadow
import scalafx.scene.input.{KeyEvent, MouseEvent}
import scalafx.scene.layout.{Pane, BorderPane}
import scalafx.scene.paint.Color
import scalafx.scene.shape.{Ellipse, Line, Polygon, Polyline, Rectangle, Shape}

abstract class TDShape
case class TDNoShape()           extends TDShape
case class Link  (shape:Line)    extends TDShape
case class Node  (shape:Ellipse) extends TDShape

object AGI extends JFXApp {

  val h = 600
  val w = 400

  var base = new MDS_Plane("mdSpace.csv","eigVals.csv")

  def projection(p:Double):Double={
    return p * 10
  }
  //
  for(i <- 1 to base.n){
    var n = new Ellipse {
      id = String.valueOf(i)
      centerX = projection(base.Xs.apply(i-1))
      centerY = projection(base.Ys.apply(i-1))
      radiusX = 10
      radiusY = 10
      stroke = Color.Black; fill = Color.Black
    }
  }

  // 平面の描画
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
        case Node(e) => e.effect = new DropShadow(10, Color.Blue); p1 = new Point2D(e.centerX(), e.centerY())
        case Link(l) => l.effect = new DropShadow(10, Color.Blue); p1 = new Point2D(l.startX(),l.startY()); p2 = new Point2D(l.endX(),l.endY())
        case TDNoShape() => ()
      }  
      p0 = new Point2D(ev.x, ev.y)
    }    

    def onDrag(ev: MouseEvent) {
      selection match {
        case Node(e) => 
          e.centerX = p1.x + ev.x - p0.x ; e.centerY = p1.y + ev.y - p0.y

        case Link(l) => 
          l.startX = p1.x + ev.x - p0.x ; l.startY = p1.y + ev.y - p0.y
          l.endX   = p2.x + ev.x - p0.x ; l.endY   = p2.y + ev.y - p0.y
        case TDNoShape() => ()
      }
    }

    def onReleased(ev: MouseEvent) {
      selection match {
        case Node(e) => base.update(ev.x,ev.y)
      }
    }

  }

  mouseAction.pressHandlers   += (("Select", SelectControl.onPress))
  mouseAction.dragHandlers    += (("Select", SelectControl.onDrag))
  mouseAction.releaseHandlers += (("Select", SelectControl.onReleased))

  val drawingPane = new Pane { }

  val shapeGroup = new ToggleGroup()

  shapeGroup.selectedToggle.onChange {
    //val id = shapeGroup.selectedToggle().asInstanceOf[javafx.scene.control.ToggleButton].id()
    //mouseAction.set(id)
  }
  mouseAction.set("Select")

  val shapeTools: List[ToggleButton] = List(
    new ToggleButton {
      id = "Select"; text = "選択"
      graphic = new Rectangle { width = 0; height = 32; fill = Color.Transparent }
      toggleGroup = shapeGroup
      minWidth = 40; minHeight = 40
    }
  )

  drawingPane.onKeyPressed = { (ev: KeyEvent) =>
    println(f"Key ${ev.code.name} was pressed")
  }

  stage = new PrimaryStage {
    title = "AGI"
    scene = new Scene(w, h) {
      root = new BorderPane {
        top = new ToolBar { content = shapeTools }
        center = drawingPane
      }
    }
  }
  drawingPane.requestFocus()
}
