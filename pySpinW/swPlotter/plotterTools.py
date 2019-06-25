import vtk

class plotterTools:

    def __init__(self, renderer):
        self.renderer = renderer
        self.MouseInteractor = plotterTools.MouseInteractor(self)

    class textBase:
        def __init__(self, pt):
            self.renderer = pt.renderer
            self.actor = None
            self.colors = vtk.vtkNamedColors()

            self.__font_size = 24
            self.__textColor = 'Tomato'
            self.__lineSpacing = 0.8
            self.__bold = False
            self.__italic = False
            self.__shadow = False

            self.__text = None
            self.__position = None
            self.__alignment = None

        @property
        def font_size(self):
            return self.__font_size

        @font_size.setter
        def font_size(self, value):
            self.__font_size = value
            self.__update()

        @property
        def textColor(self):
            return self.__textColor

        @textColor.setter
        def textColor(self, value):
            self.__textColor = value
            self.__update()

        @property
        def lineSpacing(self):
            return self.__lineSpacing

        @lineSpacing.setter
        def lineSpacing(self, value):
            self.__lineSpacing = value
            self.__update()

        @property
        def bold(self):
            return self.__bold

        @bold.setter
        def bold(self, value):
            self.__bold = value
            self.__update()

        @property
        def italic(self):
            return self.__italic

        @italic.setter
        def italic(self, value):
            self.__italic = value
            self.__update()

        @property
        def shadow(self):
            return self.__shadow

        @shadow.setter
        def shadow(self, value):
            self.__shadow = value
            self.__update()

        def addText(self, text, positionX, positionY, align='lt'):
            self.__text = text
            self.__position = [positionX, positionY]
            self.__alignment = align
            textActorL = self.textActor(text, align)
            textActorL.GetPositionCoordinate().SetValue(positionX, positionY)
            self.actor = textActorL
            self.renderer.AddActor2D(self.actor)

        def removeText(self):
            self.renderer.RemoveActor(self.actor)
            self.actor = None

        def __update(self):
            if self.__text is None:
                return
            self.removeText()
            self.addText(self.__text, self.__position[0], self.__position[1], self.__alignment)

        def textActor(self, text, align='lt'):

            def alignment(tprop, align):
                if align[0] is 'l':
                    tprop.SetJustificationToLeft()
                elif align[0] is 'r':
                    tprop.SetJustificationToRight()
                elif align[0] is 'c':
                    tprop.SetJustificationToCenter()
                else:
                    tprop.SetJustificationToLeft()
                if len(align) > 1:
                    if align[1] is 't':
                        tprop.SetVerticalJustificationToTop()
                    elif align[1] is 'b':
                        tprop.SetVerticalJustificationToBottom()
                    elif align[1] is 'c':
                        tprop.SetVerticalJustificationToCenter()
                    else:
                        tprop.SetVerticalJustificationToTop()
                else:
                    tprop.SetVerticalJustificationToTop()

            # Create the text mappers and the associated Actor2Ds.
            # The font and text properties (except justification) are the same for
            # each single line mapper. Let's create a common text property object
            singleLineTextProp = vtk.vtkTextProperty()
            singleLineTextProp.SetFontSize(self.font_size)
            singleLineTextProp.SetFontFamilyToArial()


            # The font and text properties (except justification) are the same for
            # each multi line mapper. Let's create a common text property object
            multiLineTextProp = vtk.vtkTextProperty()
            multiLineTextProp.ShallowCopy(singleLineTextProp)

            if self.bold:
                multiLineTextProp.BoldOn()
            else:
                multiLineTextProp.BoldOff()

            if self.italic:
                multiLineTextProp.ItalicOn()
            else:
                multiLineTextProp.ItalicOff()

            if self.shadow:
                multiLineTextProp.ShadowOn()
            else:
                multiLineTextProp.ShadowOff()

            multiLineTextProp.SetLineSpacing(self.lineSpacing)

            # The text is on multiple lines and left- and top-justified.
            textMapperL = vtk.vtkTextMapper()
            textMapperL.SetInput(text)
            tprop = textMapperL.GetTextProperty()
            tprop.ShallowCopy(multiLineTextProp)
            alignment(tprop, align)
            tprop.SetColor(self.colors.GetColor3d(self.textColor))

            textActorL = vtk.vtkActor2D()
            textActorL.SetMapper(textMapperL)
            textActorL.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
            return textActorL

    class tooltipPlotter(textBase):
        def __init__(self, pt):
            plotterTools.textBase.__init__(self, pt)
            self.lineSpacing = 1.05

        def setTooltip(self, text, positionX=0.05, positionY=0.95):

            if self._textBase__text is not None:
                self.removeText()
                self._textBase__text = None

            self.actor = self.textActor(text)
            self._textBase__text = text
            self.actor.GetPositionCoordinate().SetValue(positionX, positionY)
            self.renderer.AddActor2D(self.actor)

        def getTooltip(self):
            return self._textBase__text

    class MouseInteractor(vtk.vtkInteractorStyleTrackballCamera):
        def __init__(self, pt):
            self.parent = pt
            self.renderer = pt.interface.renderer
            self.interactor = pt.interface
            self.tooltipRenderer = plotterTools.tooltipPlotter(self)

        def addActors(self, actors):
            self.actorsList = actors



    class MouseClicker(MouseInteractor):

        def __init__(self, pt, type):
            plotterTools.MouseInteractor.__init__(self, pt)
            # options = ['highlightClick', 'shiftLeftClick']

            # For right click
            self.LastPickedActor = None
            self.LastPickedProperty = vtk.vtkProperty()
            self.actorsList = None

            if 'highlightClick' in type:
                self.AddObserver("RightButtonPressEvent", self.rightButtonPressEvent)

            # For left click...
            self.actor1 = None
            self.actor2 = None
            self.working = False
            self.TempActors = []

            if 'shiftLeftClick' in type:
                self.AddObserver("LeftButtonPressEvent", self.leftButtonPressEvent)

        def leftButtonPressEvent(self, obj, event):
            shitKey = self.interactor.iren.GetShiftKey()

            atomList = [atom['actor'] for atom in self.actorsList['atoms']]
            atomInds = [atom['index'] for atom in self.actorsList['atoms']]
            atomTran = [atom['unitCell'] for atom in self.actorsList['atoms']]
            indexUnitCell = [atom['indexUnitCell'] for atom in self.actorsList['atoms']]

            if shitKey and (self.actorsList is not None):

                clickPos = self.GetInteractor().GetEventPosition()

                picker = vtk.vtkPropPicker()
                picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())

                thisActor = picker.GetActor()

                if thisActor in [atom['actor'] for atom in self.actorsList['atoms']]:
                    if self.working:
                        self.actor2 = thisActor
                        self.working = False
                        self.TempActors.append(
                        self.parent.plotter.drawBond(atomInds[atomList.index(self.actor1)],
                                                     atomInds[atomList.index(self.actor2)]))

                        tempActors = self.parent.plotter.drawSymBonds(atomInds[atomList.index(self.actor1)],
                                                atomInds[atomList.index(self.actor2)],
                                                atomTran[atomList.index(self.actor2)] - atomTran[atomList.index(self.actor1)])
                        self.TempActors += tempActors

                    else:
                        self.actor1 = thisActor
                        self.working = True
                else:
                    for actor in self.TempActors:
                        self.parent.plotter.rmActor(actor)
                    self.TempActors = []
            else:
                if self.working:
                    self.working = False

            self.OnLeftButtonDown()
            return

        def rightButtonPressEvent(self, obj, event):

            if self.actorsList is None:
                return

            clickPos = self.GetInteractor().GetEventPosition()

            picker = vtk.vtkPropPicker()
            picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())

            # get the new
            self.NewPickedActor = picker.GetActor()

            # If something was selected
            if self.NewPickedActor:

                tooltipString = str('Unknown actor :-(')

                # An atom has been selected.
                if self.NewPickedActor in [atom['actor'] for atom in self.actorsList['atoms']]:
                    # If we picked something before, reset its property
                    if self.LastPickedActor:
                        self.LastPickedActor.GetProperty().DeepCopy(self.LastPickedProperty)

                    atomLabels = [atom['label'] for atom in self.actorsList['atoms']]
                    atomList   = [atom['actor'] for atom in self.actorsList['atoms']]
                    tooltipString = str(atomLabels[atomList.index(self.NewPickedActor)])

                # A lattice has been selected
                elif self.NewPickedActor in [lat['actor'] for lat in self.actorsList['lattice']]:

                    if self.LastPickedActor:
                        self.LastPickedActor.GetProperty().DeepCopy(self.LastPickedProperty)

                    latLabels = [lat['label'] for lat in self.actorsList['lattice']]
                    latList   = [lat['actor'] for lat in self.actorsList['lattice']]
                    tooltipString = str(latLabels[latList.index(self.NewPickedActor)])

                self.tooltipRenderer.setTooltip(tooltipString)
                # Save the property of the picked actor so that we can
                # restore it next time
                self.LastPickedProperty.DeepCopy(self.NewPickedActor.GetProperty())
                # Highlight the picked actor by changing its properties
                self.NewPickedActor.GetProperty().SetColor(1.0, 0.0, 0.0)
                self.NewPickedActor.GetProperty().SetDiffuse(1.0)
                self.NewPickedActor.GetProperty().SetSpecular(0.0)

                # save the last picked actor
                self.LastPickedActor = self.NewPickedActor
            else:
                if self.LastPickedActor is not None:
                    self.LastPickedActor.GetProperty().DeepCopy(self.LastPickedProperty)
                    self.LastPickedActor = None
                    if self.tooltipRenderer.getTooltip() is not None:
                        self.tooltipRenderer.removeText()
                        self.tooltipRenderer._textBase__text = None

            self.OnRightButtonDown()
            return
