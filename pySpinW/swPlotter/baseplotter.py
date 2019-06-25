import pyvista
import vtk
import numpy as np
from matplotlib import rc_params

class plotter:
    def __init__(self, spinWobj, interface, extent=[1, 1, 1]):
        self.spinwObj = spinWobj
        self.interface = interface
        self.componentsShow = {'lattice': [],
                               'atoms': [],
                               'bonds': []
                               }
        self.__extent = extent

        self.__lattice_drawn = False
        self.__lattice_redraw = False
        self.__atoms_drawn = False
        self.__atoms__redraw = False
        self.__bonds_drawn = False
        self.__bonds_redraw = False


    @property
    def extent(self):
        return self.__extent
    @extent.setter
    def extent(self, value):
        self.__extent = value
        self.update_extent()


    @property
    def latticeDrawn(self):
        if self.componentsShow['lattice']:
            self.__lattice_drawn = True
        else:
            self.__lattice_drawn = False
        return self.__lattice_drawn


    @property
    def atomsDrawn(self):
        if self.componentsShow['atoms']:
            self.__atoms_drawn = True
        else:
            self.__atoms_drawn = False
        return self.__atoms_drawn

    @property
    def bondsDrawn(self):
        if self.componentsShow['bonds']:
            self.__bonds_drawn = True
        else:
            self.__bonds_drawn = False
        return self.__bonds_drawn


    def frameCycler(self):

            #Do we have a lattice?
            if self.latticeDrawn:
                self.__lattice_redraw = True
                self.rm_lattice()

            #Do we have atoms?
            if self.atomsDrawn:
                self.rm_atoms()
                self.__atoms__redraw = True

            # Do we have bonds?
            if self.bondsDrawn:
                self.rm_bonds()
                self.__bonds_redraw = True
            if self.__lattice_redraw:
                self.add_lattice()
                self.__lattice_redraw = False
            if self.__atoms__redraw:
                self.add_atoms()
                self.__atoms__redraw = False
            if self.bondsDrawn:
                self.add_bonds()
                self.__bonds_redraw =False


    def update_extent(self, extent=None):
        if extent is not None:
            self.__extent = extent
        self.frameCycler()

    def add_lattice(self):
        """ add a lattice to the pyqt frame """
        if not self.componentsShow['lattice']:
            # bv = np.linalg.inv(self.spinwObj.crystal.get_reciprocal_cell())
            bv = self.spinwObj.basisvector()

            extent = self.extent
            lines = []
            tempStore = []

            def orderPoints(ptA, ptB):
                points = np.mod(np.array([ptA, ptB]), 2)
                idx = np.argsort(np.sum(points * [3, 2, 1], axis=1))
                points = points[idx,:]
                return (points[0,:], points[1,:])



            # Add major lattice points
            for z in np.arange(extent[2]):
                for y in np.arange(extent[1]):
                    for x in np.arange(extent[0]):
                        lines.append(pyvista.Line(pointa=np.matmul(bv, [x, y, z] + 1E-5*(np.random.rand(3)-0.5)), pointb=np.matmul(bv, [x + 1, y, z] + 1E-5*(np.random.rand(3)-0.5))))
                        startPt, endPt = orderPoints([x, y, z], [x+1, y, z])
                        tempStore.append({'actor': None, 'label': 'Lattice\nCell: [%i, %i, %i]\nVector: [%i, %i, %i] - [%i, %i, %i]' % (x+1, y+1, z+1, startPt[0], startPt[1], startPt[2], endPt[0], endPt[1], endPt[2])})
                        lines.append(pyvista.Line(pointa=np.matmul(bv, [x, y, z] + 1E-5*(np.random.rand(3)-0.5)), pointb=np.matmul(bv, [x, y + 1, z] + 1E-5*(np.random.rand(3)-0.5))))
                        startPt, endPt = orderPoints([x, y, z], [x, y+1, z])
                        tempStore.append({'actor': None, 'label': 'Lattice\nCell: [%i, %i, %i]\nVector: [%i, %i, %i] - [%i, %i, %i]' % (x+1, y+1, z+1, startPt[0], startPt[1], startPt[2], endPt[0], endPt[1], endPt[2])})
                        lines.append(pyvista.Line(pointa=np.matmul(bv, [x, y, z] + 1E-5*(np.random.rand(3)-0.5)), pointb=np.matmul(bv, [x, y, z + 1] + 1E-5*(np.random.rand(3)-0.5))))
                        startPt, endPt = orderPoints([x, y, z], [x, y, z+1])
                        tempStore.append({'actor': None, 'label': 'Lattice\nCell: [%i, %i, %i]\nVector: [%i, %i, %i] - [%i, %i, %i]' % (x+1, y+1, z+1, startPt[0], startPt[1], startPt[2], endPt[0], endPt[1], endPt[2])})

                        lines.append(pyvista.Line(pointa=np.matmul(bv, [x + 1, y + 1, z + 1] + 1E-5*(np.random.rand(3)-0.5)),
                                               pointb=np.matmul(bv, [x + 1, y, z + 1] + 1E-5*(np.random.rand(3)-0.5))))
                        startPt, endPt = orderPoints([x+1, y+1, z+1], [x + 1, y, z+1])
                        tempStore.append({'actor': None, 'label': 'Lattice\nCell: [%i, %i, %i]\nVector: [%i, %i, %i] - [%i, %i, %i]' % (x+1, y+1, z+1, startPt[0], startPt[1], startPt[2], endPt[0], endPt[1], endPt[2])})
                        lines.append(pyvista.Line(pointa=np.matmul(bv, [x + 1, y + 1, z + 1] + 1E-5*(np.random.rand(3)-0.5)),
                                               pointb=np.matmul(bv, [x + 1, y + 1, z] + 1E-5*(np.random.rand(3)-0.5))))
                        startPt, endPt = orderPoints([x + 1, y + 1, z + 1], [x + 1, y+1, z])
                        tempStore.append({'actor': None,
                                          'label': 'Lattice\nCell: [%i, %i, %i]\nVector: [%i, %i, %i] - [%i, %i, %i]' % (x+1, y+1, z+1, startPt[0], startPt[1], startPt[2], endPt[0], endPt[1], endPt[2])})
                        lines.append(pyvista.Line(pointa=np.matmul(bv, [x + 1, y + 1, z + 1] + 1E-5*(np.random.rand(3)-0.5)),
                                               pointb=np.matmul(bv, [x, y + 1, z + 1] + 1E-5*(np.random.rand(3)-0.5))))
                        startPt, endPt = orderPoints([x + 1, y + 1, z + 1], [x, y+1, z + 1])
                        tempStore.append({'actor': None, 'label': 'Lattice\nCell: [%i, %i, %i]\nVector: [%i, %i, %i] - [%i, %i, %i]' % (x+1, y+1, z+1, startPt[0], startPt[1], startPt[2], endPt[0], endPt[1], endPt[2])})
            # Draw x lines
            for x in np.arange(extent[0]):
                y = extent[1]
                z = 0
                lines.append(pyvista.Line(pointa=np.matmul(bv, [x, y, z] + 1E-5*(np.random.rand(3)-0.5)), pointb=np.matmul(bv, [x + 1, y, z] + 1E-5*(np.random.rand(3)-0.5))))
                startPt, endPt = orderPoints([x, y, z], [x + 1, y, z])
                tempStore.append({'actor': None, 'label': 'Lattice\nCell: [%i, %i, %i]\nVector: [%i, %i, %i] - [%i, %i, %i]' %
                                                          (x+1, y, z+1, startPt[0], startPt[1], startPt[2], endPt[0], endPt[1], endPt[2])})
                z = extent[2]
                y = 0
                lines.append(pyvista.Line(pointa=np.matmul(bv, [x, y, z] + 1E-5*(np.random.rand(3)-0.5)), pointb=np.matmul(bv, [x + 1, y, z] + 1E-5*(np.random.rand(3)-0.5))))
                startPt, endPt = orderPoints([x, y, z], [x + 1, y, z])
                tempStore.append({'actor': None, 'label': 'Lattice\nCell: [%i, %i, %i]\nVector: [%i, %i, %i] - [%i, %i, %i]' %
                                                          (x+1, y+1, z, startPt[0], startPt[1], startPt[2], endPt[0], endPt[1], endPt[2])})

            # Draw y lines
            for y in np.arange(extent[1]):
                x = 0
                z = extent[2]
                lines.append(pyvista.Line(pointa=np.matmul(bv, [x, y, z] + 1E-5*(np.random.rand(3)-0.5)), pointb=np.matmul(bv, [x, y + 1, z] + 1E-5*(np.random.rand(3)-0.5))))
                startPt, endPt = orderPoints([x, y, z], [x, y + 1, z])
                tempStore.append({'actor': None, 'label': 'Lattice\nCell: [%i, %i, %i]\nVector: [%i, %i, %i] - [%i, %i, %i]' %
                                                          (x+1, y+1, z, startPt[0], startPt[1], startPt[2], endPt[0], endPt[1], endPt[2])})
                x = extent[0]
                z = 0
                lines.append(pyvista.Line(pointa=np.matmul(bv, [x, y, z] + 1E-5*(np.random.rand(3)-0.5)), pointb=np.matmul(bv, [x, y + 1, z] + 1E-5*(np.random.rand(3)-0.5))))
                startPt, endPt = orderPoints([x, y, z], [x, y + 1, z])
                tempStore.append({'actor': None, 'label': 'Lattice\nCell: [%i, %i, %i]\nVector: [%i, %i, %i] - [%i, %i, %i]' %
                                                          (x, y+1, z+1, startPt[0], startPt[1], startPt[2], endPt[0], endPt[1], endPt[2])})

            # Draw y lines
            for z in np.arange(extent[2]):
                x = 0
                y = extent[1]
                lines.append(pyvista.Line(pointa=np.matmul(bv, [x, y, z] + 1E-5*(np.random.rand(3)-0.5)), pointb=np.matmul(bv, [x, y, z + 1] + 1E-5*(np.random.rand(3)-0.5))))
                startPt, endPt = orderPoints([x, y, z], [x, y, z + 1])
                tempStore.append({'actor': None, 'label': 'Lattice\nCell: [%i, %i, %i]\nVector: [%i, %i, %i] - [%i, %i, %i]' %
                                                          (x+1, y, z+1, startPt[0], startPt[1], startPt[2], endPt[0], endPt[1], endPt[2])})
                x = extent[0]
                y = 0
                lines.append(pyvista.Line(pointa=np.matmul(bv, [x, y, z] + 1E-5*(np.random.rand(3)-0.5)), pointb=np.matmul(bv, [x, y, z + 1] + 1E-5*(np.random.rand(3)-0.5))))
                startPt, endPt = orderPoints([x, y, z], [x, y, z + 1])
                tempStore.append({'actor': None, 'label': 'Lattice\nCell: [%i, %i, %i]\nVector: [%i, %i, %i] - [%i, %i, %i]' %
                                                          (x, y+1, z+1, startPt[0], startPt[1], startPt[2], endPt[0], endPt[1], endPt[2])})

            for n, line in enumerate(lines):
                tempStore[n]['actor'] = self.interface.add_mesh(line, color='black')
                self.componentsShow['lattice'].append(tempStore[n].copy())
            self.interface.reset_camera()

    def rm_lattice(self):

        # Remove the mesh from the plot
        for line in self.componentsShow['lattice']:
            self.rmActor(line['actor'])
        # Remove the meshes from the list
        self.componentsShow['lattice'] = []

    def add_atoms(self):
        """ add a sphere to the pyqt frame """
        atom = self.spinwObj.atom()

        extent = np.array(self.extent)

        rangelu = np.ceil(extent)
        cells = np.flip(np.flip(np.mgrid[0:rangelu[2]+1:1, 0:rangelu[1]+1:1, 0:rangelu[0]+1:1].reshape((3, -1))), axis=1)
        pos = []
        for column in cells.T:
            pos.append(np.add(atom['r'], column.reshape(3, 1)))
        pos = np.array(pos).reshape((3, -1))
        nCell = np.prod(rangelu+1)
        inRange = np.all((pos > [[0], [0], [0]]-10*np.finfo(float).eps) &
                         (pos < extent.reshape((3, 1))+10*np.finfo(float).eps), axis=0)
        atomIDX = np.repeat(atom['idx'], nCell)

        pos = pos[:, inRange]
        atomIDX = atomIDX[inRange]

        t = rc_params()['axes.prop_cycle']
        col = [thisCol for idx, thisCol in zip(np.unique(atom['idx']), t)]

        if self.latticeDrawn:
            for line in self.componentsShow['lattice']:
                line['actor'].SetVisibility(0)

        for ind, p in enumerate(pos.T):
            sphere = pyvista.Sphere(radius=0.2, center=np.matmul(self.spinwObj.basisvector(), p))
            mesh = self.interface.add_mesh(sphere, color=col[int(atomIDX[ind])-1]['color'])
            self.componentsShow['atoms'].append({'actor': mesh,
                                                 'label': 'Atom\nLabel:{:s}\nAtomic position:[{:.3f} {:.3f} {:.3f}]\nUnit cell:[{:d} {:d} {:d}]'.format(
                                                     self.spinwObj.unit_cell['label'][int(atomIDX[ind])-1],
                                                     p[0] - cells[0, int(np.mod(ind, nCell))],
                                                     p[1] - cells[2, int(np.mod(ind, nCell))],
                                                     p[2] - cells[2, int(np.mod(ind, nCell))],
                                                     int(cells[0, int(np.mod(ind, nCell))]),
                                                     int(cells[1, int(np.mod(ind, nCell))]),
                                                     int(cells[2, int(np.mod(ind, nCell))])),
                                                 'position': pos,
                                                 'positionUnitCell': p - cells[:, int(np.mod(ind, nCell))],
                                                 'unitCell': cells[:, int(np.mod(ind, nCell))],
                                                 'symbol': self.spinwObj.unit_cell['label'][int(atomIDX[ind])-1],
                                                 'indexUnitCell': atomIDX[ind],
                                                 'index': ind})

        if not self.componentsShow['lattice']:
            redoLattice = False
        else:
            self.rm_lattice()
            redoLattice = True

        if self.latticeDrawn:
            for line in self.componentsShow['lattice']:
                line['actor'].SetVisibility(1)

        if redoLattice:
            self.add_lattice()

    def rm_atoms(self):

        # Remove the mesh from the plot
        for atom in self.componentsShow['atoms']:
            self.rmActor(atom['actor'])
        # Remove the meshes from the list
        self.componentsShow['atoms'] = []


    def rmActor(self, actor):
        self.interface.remove_actor(actor)

    def add_bonds(self):

        # bv = np.linalg.inv(self.spinwObj.crystal.get_reciprocal_cell())

        bv = self.spinwObj.basisvector()

        # if not self.componentsShow['lattice']:
        #     redoLattice = False
        # else:
        #     self.rm_lattice()
        #     redoLattice = True

        bond_idx = 0

        start = self.spinwObj.crystal[self.spinwObj.bonds['atom1'].astype(int)][self.spinwObj.bonds['idx']==bond_idx].positions
        end = self.spinwObj.crystal[self.spinwObj.bonds['atom2'].astype(int)][self.spinwObj.bonds['idx']==bond_idx].positions
        dl = np.matmul(bv,self.spinwObj.bonds['dl'][self.spinwObj.bonds['idx']==bond_idx].reshape((-1,3,1))).reshape((-1, 3))


        centres = (start + end + dl) / 2
        direction = end - start - dl
        # direction = np.matmul(direction, bv.transpose())
        # lengths = np.sqrt(np.sum(np.matmul(bv, direction.transpose())**2, axis=0))
        lengths = np.sqrt(np.sum(direction.transpose() ** 2, axis=0))

        bonds = []

        for ind, point in enumerate(centres):
            cylinder = pyvista.Cylinder(center=point+1E-5*(np.random.rand()-0.5), direction=direction[ind, :], height=lengths[ind], radius=0.025)
            # if self.spinwObj.bonds['type'][ind] == 0:
            bonds.append(self.interface.add_mesh(cylinder, color='red'))
            # elif self.spinwObj.bonds['type'][ind] == 1:
            #     bonds.append(self.interface.add_mesh(cylinder, color='green'))
            self.componentsShow['bonds'].append(cylinder)

        self.componentsShow['bonds'] = bonds
        # if redoLattice:
        #     self.add_lattice()


    def drawBond(self, atom1, atom2, dl=[0, 0, 0], color='red'):

        bv = np.linalg.inv(self.spinwObj.crystal.get_reciprocal_cell())
        extent = np.array(self.extent)
        atomsL = self.spinwObj.crystal.repeat((extent[0] + 1, extent[1] + 1, extent[2] + 1))
        start = atomsL[int(atom1)].position
        end = atomsL[int(atom2)].position

        a = np.any(np.matmul(self.spinwObj.crystal.get_reciprocal_cell(), end) + dl > extent)
        b = np.any(np.matmul(self.spinwObj.crystal.get_reciprocal_cell(), start) + dl < [0, 0, 0])
        #
        if a or b:
            return

        # dl = np.matmul(bv,self.spinwObj.bonds['dl'][self.spinwObj.bonds['idx']==bond_idx].reshape((-1,3,1))).reshape((-1, 3))
        dl = np.matmul(bv, dl)

        centres = (start + end + dl) / 2
        direction = end + dl - start

        # direction = np.matmul(direction, bv.transpose())
        # lengths = np.sqrt(np.sum(np.matmul(bv, direction.transpose())**2, axis=0))
        lengths = np.sqrt(np.sum(direction.transpose() ** 2, axis=0))

        cylinder = pyvista.Cylinder(center=centres, direction=direction, height=lengths, radius=0.025)
        actor = self.interface.add_mesh(cylinder, color=color)
        # self.frameCycler()
        return actor

    def drawSymBonds(self, atom1, atom2, dl):

        extent = np.array(self.extent)
        atomsL = self.spinwObj.crystal.repeat((extent[0] + 1, extent[1] + 1, extent[2] + 1))
        spos = atomsL.get_scaled_positions() * (extent + 1) / extent

        atom1 = np.where(np.all(spos == spos[atom1] - np.floor(spos[atom1] - dl), axis=1))[0][0]
        atom2 = np.where(np.all(spos == spos[atom2] - np.floor(spos[atom2] - dl), axis=1))[0][0]

        r = np.array([dl[0], dl[1], dl[2], atom1, atom2])
        oneBonds, newBond = self.spinwObj.crystal.info['spacegroup'].bonds(self.spinwObj.crystal.get_scaled_positions(),
                                                              self.spinwObj.basisvector(),
                                                              r)
        bonds = np.tile(oneBonds, [4, 1])
        bonds[1*len(oneBonds):2 * len(oneBonds), 0] = bonds[1*len(oneBonds):2 * len(oneBonds), 0] + 1
        bonds[2*len(oneBonds):3 * len(oneBonds), 1] = bonds[2*len(oneBonds):3 * len(oneBonds), 0] + 1
        bonds[3*len(oneBonds):4 * len(oneBonds), 2] = bonds[3*len(oneBonds):4 * len(oneBonds), 0] + 1

        actors = []

        for bond in bonds:
            if np.all(bond == r):
                #Self bond. Already plotted......
                continue
            dl = bond[0:3]
            atom1 = bond[3]
            atom2 = bond[4]
            actors.append(self.drawBond(atom1, atom2, dl, color='gold'))

        # for bond in bonds[np.all(bonds[:,0:3]==[0, 0, 0], axis=1) & newBond, :]:
        #     if np.all(bond == r):
        #         #Self bond
        #         continue
        #     dl = bond[0:3]
        #     atom1 = bond[3]
        #     atom2 = bond[4]
        #     actors.append(self.drawBond(atom1, atom2, dl, color='gold'))



        # for ex0 in np.arange(self.extent[0]+1):
        #     for ex1 in np.arange(self.extent[1]+1):
        #         for ex2 in np.arange(self.extent[2]+1):
        #             for bond in bonds[np.all(bonds[:, 0:3] == [ex0, ex1, ex2], axis=1) & newBond, :]:
        #                 if np.all(bond == r):
        #                     #Self bond
        #                     continue
        #                 dl = bond[0:3]
        #                 atom1 = bond[3]
        #                 atom2 = bond[4]
        #                 actors.append(self.drawBond(atom1, atom2, dl, color='gold'))

        return actors

    def rm_bonds(self):

        # Remove the mesh from the plot
        for mesh in self.componentsShow['bonds']:
            self.rmActor(mesh)
        # Remove the meshes from the list
        self.componentsShow['bonds'] = []

    def rotate(self, angle1=0, angle2=0, axis=None):
        # See https://gist.github.com/pangyuteng/facd430d0d9761fc67fff4ff2e5fffc3

        bv = self.spinwObj.basisvector()
        center = np.matmul(np.array([0.5, 0.5, 0.5])*self.extent, bv.transpose())

        focal_point = self.interface.camera.GetFocalPoint()
        view_up = self.interface.camera.GetViewUp()
        position = self.interface.camera.GetPosition()

        if axis is None:
            axis = [0, 0, 0]
            axis[0] = -1 * self.interface.camera.GetViewTransformMatrix().GetElement(0, 0)
            axis[1] = -1 * self.interface.camera.GetViewTransformMatrix().GetElement(0, 1)
            axis[2] = -1 * self.interface.camera.GetViewTransformMatrix().GetElement(0, 2)


        transform = vtk.vtkTransform()
        transform.Identity()

        transform.Translate(*center)
        transform.RotateWXYZ(angle1, view_up)
        transform.RotateWXYZ(angle2, axis)
        transform.Translate(*[-1 * x for x in center])

        new_position = [0, 0, 0]
        new_focal_point = [0, 0, 0]
        transform.TransformPoint(position, new_position)
        transform.TransformPoint(focal_point, new_focal_point)

        self.interface.camera.SetPosition(new_position)
        self.interface.camera.SetFocalPoint(new_focal_point)
        self.interface.camera.OrthogonalizeViewUp()
        self.interface.renderer.ResetCameraClippingRange()

    def along_a(self):
        self.interface.reset_camera()
        bv = self.spinwObj.basisvector()
        center = np.matmul(np.array([0.5, 0.5, 0.5]) * self.extent, bv.transpose())
        self.interface.reset_camera()
        self.interface.camera.SetPosition([1, 0, 0])
        self.interface.camera.SetFocalPoint(center)
        self.interface.camera.OrthogonalizeViewUp()
        self.interface.camera.Roll(0)
        self.interface.renderer.ResetCameraClippingRange()

    def along_b(self):
        self.interface.reset_camera()
        bv = self.spinwObj.basisvector()
        center = np.matmul(np.array([0.5, 0.5, 0.5]) * self.extent, bv.transpose())
        self.interface.camera.SetPosition([0, 1, 0])
        self.interface.camera.SetFocalPoint(center)
        self.interface.camera.Roll(0)
        self.interface.camera.OrthogonalizeViewUp()
        self.interface.renderer.ResetCameraClippingRange()

    def along_c(self):
        bv = self.spinwObj.basisvector()
        center = np.matmul(np.array([0.5, 0.5, 0.5]) * self.extent, bv.transpose())
        self.interface.camera.SetPosition([0, 0, 1])
        self.interface.camera.SetFocalPoint(center)
        self.interface.camera.Roll(0)
        self.interface.camera.OrthogonalizeViewUp()
        self.interface.renderer.ResetCameraClippingRange()
        self.interface.reset_camera()

        # self.interface.reset_camera()
        # bv = self.spinwObj.basisvector()
        # center = np.matmul(np.array([0.5, 0.5, 0.5]) * self.extent, bv.transpose())
        # position = self.interface.camera.GetPosition()
        # rot = self.interface.ren_win.
        # focal_point = center
        # view_up = [0, 1, 0]
        # axis = [0, 0, -1]
        # transform = vtk.vtkTransform()
        # transform.Identity()
        #
        # transform.Translate(*center)
        # transform.RotateWXYZ(0, view_up)
        # transform.RotateWXYZ(0, axis)
        # transform.Translate(*[-1 * x for x in center])
        #
        # new_position = [0, 0, 0]
        # new_focal_point = [0, 0, 0]
        # transform.TransformPoint(position, new_position)
        # transform.TransformPoint(focal_point, new_focal_point)
        #
        # self.interface.camera.SetPosition(new_position)
        # self.interface.camera.SetFocalPoint(new_focal_point)
        # self.interface.camera.OrthogonalizeViewUp()
        # self.interface.renderer.ResetCameraClippingRange()