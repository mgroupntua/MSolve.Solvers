//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;
//using MGroup.MSolve.Discretization.Mesh;
//using MGroup.Solvers.DDM.Mesh;
//using MGroup.Solvers.DDM.Output.VTK;
//using Xunit;

//namespace MGroup.Solvers.DDM.Tests.Mesh
//{
//    public static class UniformMesh3DTests
//    {
//        //[Fact]
//        public static void PlotMesh()
//        {
//            var writer = new VtkMeshWriter();

//            string path = @"C:\Users\Serafeim\Desktop\PFETIDP\meshes\mesh3D.vtk";
//            double[] minCoords = { 0, 0, 0 };
//            double[] maxCoords = { 60, 60, 60 };
//            int[] numElements = { 2, 3, 4 };
//            //var mesh = new UniformMesh3D.Builder(minCoords, maxCoords, numElements).SetMajorMinorAxis(2, 0).BuildMesh();
//            var mesh = new UniformMesh3D.Builder(minCoords, maxCoords, numElements).BuildMesh();
//            writer.WriteMesh(path, mesh, 3);
//        }

//        [Fact]
//        public static void TestEnumerateElements()
//        {
//            var mockMesh = new MockMesh2x3x4();
//            var mesh = new UniformMesh3D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
//                .BuildMesh();

//            // Check that nodes are returned in ascending order
//            int[] elemIDs = mesh.EnumerateElements().Select(pair => pair.elementID).ToArray();
//            Assert.Equal(0, elemIDs[0]);
//            for (int i = 0; i < elemIDs.Length - 1; ++i)
//            {
//                Assert.Equal(elemIDs[i] + 1, elemIDs[i + 1]);
//            }

//            // Check each returned element
//            foreach ((int elementID, int[] nodeIDs) in mesh.EnumerateElements())
//            {
//                int[] nodeIDsExpected = mockMesh.GetElementConnectivity(mockMesh.GetElementIdx(elementID));
//                for (int d = 0; d < nodeIDsExpected.Length; ++d)
//                {
//                    Assert.Equal(nodeIDsExpected[d], nodeIDs[d]);
//                }
//            }
//        }

//        [Fact]
//        public static void TestEnumerateNodes()
//        {
//            var mockMesh = new MockMesh2x3x4();
//            var mesh = new UniformMesh3D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
//                .BuildMesh();

//            // Check that nodes are returned in ascending order
//            int[] nodeIDs = mesh.EnumerateNodes().Select(pair => pair.nodeID).ToArray();
//            Assert.Equal(0, nodeIDs[0]);
//            for (int i = 0; i < nodeIDs.Length - 1; ++i)
//            {
//                Assert.Equal(nodeIDs[i] + 1, nodeIDs[i + 1]);
//            }

//            // Check each returned node
//            foreach ((int nodeID, double[] coords) in mesh.EnumerateNodes())
//            {
//                double[] coordsExpected = mockMesh.GetNodeCoordinates(mockMesh.GetNodeIdx(nodeID));
//                for (int d = 0; d < 3; ++d)
//                {
//                    Assert.Equal(coordsExpected[d], coords[d]);
//                }
//            }
//        }

//        [Fact]
//        public static void TestGetNodeID()
//        {
//            var mockMesh = new MockMesh2x3x4();
//            var mesh = new UniformMesh3D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
//                .BuildMesh();

//            for (int i = 0; i < mesh.NumNodes[0]; ++i)
//            {
//                for (int j = 0; j < mesh.NumNodes[1]; ++j)
//                {
//                    for (int k = 0; k < mesh.NumNodes[2]; ++k)
//                    {
//                        var idx = new int[] { i, j, k };
//                        int expected = mockMesh.GetNodeID(idx);
//                        int computed = mesh.GetNodeID(idx);
//                        Assert.Equal(expected, computed);
//                    }
//                }
//            }
//        }

//        [Fact]
//        public static void TestGetNodeIdx()
//        {
//            var mockMesh = new MockMesh2x3x4();
//            var mesh = new UniformMesh3D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
//                .BuildMesh();

//            for (int id = 0; id < mesh.NumNodes[0] * mesh.NumNodes[1] * mesh.NumNodes[2]; ++id)
//            {
//                int[] expected = mockMesh.GetNodeIdx(id);
//                int[] computed = mesh.GetNodeIdx(id);
//                for (int d = 0; d < 3; ++d)
//                {
//                    Assert.Equal(expected[d], computed[d]);
//                }
//            }
//        }

//        [Fact]
//        public static void TestGetNodeCoordinates()
//        {
//            var mockMesh = new MockMesh2x3x4();
//            var mesh = new UniformMesh3D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
//                .BuildMesh();

//            for (int i = 0; i < mesh.NumNodes[0]; ++i)
//            {
//                for (int j = 0; j < mesh.NumNodes[1]; ++j)
//                {
//                    for (int k = 0; k < mesh.NumNodes[2]; ++k)
//                    {
//                        var idx = new int[] { i, j, k };
//                        double[] expected = mockMesh.GetNodeCoordinates(idx);
//                        double[] computed = mesh.GetNodeCoordinates(idx);
//                        for (int d = 0; d < 3; ++d)
//                        {
//                            Assert.Equal(expected[d], computed[d]);
//                        }
//                    }
//                }
//            }
//        }

//        [Fact]
//        public static void TestGetElementID()
//        {
//            var mockMesh = new MockMesh2x3x4();
//            var mesh = new UniformMesh3D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
//                .BuildMesh();

//            for (int i = 0; i < mesh.NumElements[0]; ++i)
//            {
//                for (int j = 0; j < mesh.NumElements[1]; ++j)
//                {
//                    for (int k = 0; k < mesh.NumElements[2]; ++k)
//                    {
//                        var idx = new int[] { i, j, k };
//                        int expected = mockMesh.GetElementID(idx);
//                        int computed = mesh.GetElementID(idx);
//                        Assert.Equal(expected, computed);
//                    }
//                }
//            }
//        }

//        [Fact]
//        public static void TestGetElementIdx()
//        {
//            var mockMesh = new MockMesh2x3x4();
//            var mesh = new UniformMesh3D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
//                .BuildMesh();

//            for (int id = 0; id < mesh.NumElements[0] * mesh.NumElements[1] * mesh.NumElements[2]; ++id)
//            {
//                int[] expected = mockMesh.GetElementIdx(id);
//                int[] computed = mesh.GetElementIdx(id);
//                for (int d = 0; d < 3; ++d)
//                {
//                    Assert.Equal(expected[d], computed[d]);
//                }
//            }
//        }

//        [Fact]
//        public static void TestGetElementConductivity()
//        {
//            var mockMesh = new MockMesh2x3x4();
//            var mesh = new UniformMesh3D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
//                .BuildMesh();

//            for (int k = 0; k < mesh.NumElements[2]; ++k)
//            {
//                for (int j = 0; j < mesh.NumElements[1]; ++j)
//                {
//                    for (int i = 0; i < mesh.NumElements[0]; ++i)
//                    {
//                        var idx = new int[] { i, j, k };
//                        int[] expected = mockMesh.GetElementConnectivity(idx);
//                        int[] computed = mesh.GetElementConnectivity(idx);
//                        for (int d = 0; d < 8; ++d)
//                        {
//                            Assert.Equal(expected[d], computed[d]);
//                        }
//                    }
//                }
//            }
//        }

//        private class MockMesh2x3x4 : IStructuredMesh
//        {
//            public CellType CellType => CellType.Hexa8;

//            public double[] MinCoordinates => new double[] { 0.0, 0.0, 0.0 };

//            public double[] MaxCoordinates => new double[] { 2.0, 3.0, 4.0 };

//            public int[] NumElements => new int[] { 2, 3, 4 };

//            public int[] NumNodes => new int[] { 3, 4, 5 };

//            public int NumElementsTotal => 24;

//            public int NumNodesTotal => 60;

//            public int NumNodesPerElement => 8;

//            public IEnumerable<(int elementID, int[] nodeIDs)> EnumerateElements()
//            {
//                throw new NotImplementedException();
//            }

//            public IEnumerable<(int nodeID, double[] coordinates)> EnumerateNodes()
//            {
//                throw new NotImplementedException();
//            }

//            public int[] GetElementConnectivity(int[] elementIdx)
//            {
//                var elementNodes = new int[2, 3, 4][];
//                elementNodes[0, 0, 0] = new int[] { 0, 1, 4, 3, 12, 13, 16, 15 };
//                elementNodes[1, 0, 0] = new int[] { 1, 2, 5, 4, 13, 14, 17, 16 };
//                elementNodes[0, 1, 0] = new int[] { 3, 4, 7, 6, 15, 16, 19, 18 };
//                elementNodes[1, 1, 0] = new int[] { 4, 5, 8, 7, 16, 17, 20, 19 };
//                elementNodes[0, 2, 0] = new int[] { 6, 7, 10, 9, 18, 19, 22, 21 };
//                elementNodes[1, 2, 0] = new int[] { 7, 8, 11, 10, 19, 20, 23, 22 };

//                elementNodes[0, 0, 1] = new int[] { 12, 13, 16, 15, 24, 25, 28, 27 };
//                elementNodes[1, 0, 1] = new int[] { 13, 14, 17, 16, 25, 26, 29, 28 };
//                elementNodes[0, 1, 1] = new int[] { 15, 16, 19, 18, 27, 28, 31, 30 };
//                elementNodes[1, 1, 1] = new int[] { 16, 17, 20, 19, 28, 29, 32, 31 };
//                elementNodes[0, 2, 1] = new int[] { 18, 19, 22, 21, 30, 31, 34, 33 };
//                elementNodes[1, 2, 1] = new int[] { 19, 20, 23, 22, 31, 32, 35, 34 };

//                elementNodes[0, 0, 2] = new int[] { 24, 25, 28, 27, 36, 37, 40, 39 };
//                elementNodes[1, 0, 2] = new int[] { 25, 26, 29, 28, 37, 38, 41, 40 };
//                elementNodes[0, 1, 2] = new int[] { 27, 28, 31, 30, 39, 40, 43, 42 };
//                elementNodes[1, 1, 2] = new int[] { 28, 29, 32, 31, 40, 41, 44, 43 };
//                elementNodes[0, 2, 2] = new int[] { 30, 31, 34, 33, 42, 43, 46, 45 };
//                elementNodes[1, 2, 2] = new int[] { 31, 32, 35, 34, 43, 44, 47, 46 };

//                elementNodes[0, 0, 3] = new int[] { 36, 37, 40, 39, 48, 49, 52, 51 };
//                elementNodes[1, 0, 3] = new int[] { 37, 38, 41, 40, 49, 50, 53, 52 };
//                elementNodes[0, 1, 3] = new int[] { 39, 40, 43, 42, 51, 52, 55, 54 };
//                elementNodes[1, 1, 3] = new int[] { 40, 41, 44, 43, 52, 53, 56, 55 };
//                elementNodes[0, 2, 3] = new int[] { 42, 43, 46, 45, 54, 55, 58, 57 };
//                elementNodes[1, 2, 3] = new int[] { 43, 44, 47, 46, 55, 56, 59, 58 };

//                return elementNodes[elementIdx[0], elementIdx[1], elementIdx[2]];
//            }

//            public int GetElementID(int[] elementIdx)
//            {
//                var elemIDs = new int[2, 3, 4];
//                elemIDs[0, 0, 0] = 0;
//                elemIDs[1, 0, 0] = 1;
//                elemIDs[0, 1, 0] = 2;
//                elemIDs[1, 1, 0] = 3;
//                elemIDs[0, 2, 0] = 4;
//                elemIDs[1, 2, 0] = 5;

//                elemIDs[0, 0, 1] = 6;
//                elemIDs[1, 0, 1] = 7;
//                elemIDs[0, 1, 1] = 8;
//                elemIDs[1, 1, 1] = 9;
//                elemIDs[0, 2, 1] = 10;
//                elemIDs[1, 2, 1] = 11;

//                elemIDs[0, 0, 2] = 12;
//                elemIDs[1, 0, 2] = 13;
//                elemIDs[0, 1, 2] = 14;
//                elemIDs[1, 1, 2] = 15;
//                elemIDs[0, 2, 2] = 16;
//                elemIDs[1, 2, 2] = 17;

//                elemIDs[0, 0, 3] = 18;
//                elemIDs[1, 0, 3] = 19;
//                elemIDs[0, 1, 3] = 20;
//                elemIDs[1, 1, 3] = 21;
//                elemIDs[0, 2, 3] = 22;
//                elemIDs[1, 2, 3] = 23;

//                return elemIDs[elementIdx[0], elementIdx[1], elementIdx[2]];
//            }

//            public int[] GetElementIdx(int elementID)
//            {
//                var elemIndices = new int[,]
//                {
//                    { 0, 0, 0 },
//                    { 1, 0, 0 },
//                    { 0, 1, 0 },
//                    { 1, 1, 0 },
//                    { 0, 2, 0 },
//                    { 1, 2, 0 },

//                    { 0, 0, 1 },
//                    { 1, 0, 1 },
//                    { 0, 1, 1 },
//                    { 1, 1, 1 },
//                    { 0, 2, 1 },
//                    { 1, 2, 1 },

//                    { 0, 0, 2 },
//                    { 1, 0, 2 },
//                    { 0, 1, 2 },
//                    { 1, 1, 2 },
//                    { 0, 2, 2 },
//                    { 1, 2, 2 },

//                    { 0, 0, 3 },
//                    { 1, 0, 3 },
//                    { 0, 1, 3 },
//                    { 1, 1, 3 },
//                    { 0, 2, 3 },
//                    { 1, 2, 3 },
//                };

//                return new int[] { elemIndices[elementID, 0], elemIndices[elementID, 1], elemIndices[elementID, 2] };
//            }

//            public double[] GetNodeCoordinates(int[] nodeIdx)
//            {
//                var nodeCoords = new double[3, 4, 5][];
//                nodeCoords[0, 0, 0] = new double[] { 0, 0, 0 };
//                nodeCoords[1, 0, 0] = new double[] { 1, 0, 0 };
//                nodeCoords[2, 0, 0] = new double[] { 2, 0, 0 };
//                nodeCoords[0, 1, 0] = new double[] { 0, 1, 0 };
//                nodeCoords[1, 1, 0] = new double[] { 1, 1, 0 };
//                nodeCoords[2, 1, 0] = new double[] { 2, 1, 0 };
//                nodeCoords[0, 2, 0] = new double[] { 0, 2, 0 };
//                nodeCoords[1, 2, 0] = new double[] { 1, 2, 0 };
//                nodeCoords[2, 2, 0] = new double[] { 2, 2, 0 };
//                nodeCoords[0, 3, 0] = new double[] { 0, 3, 0 };
//                nodeCoords[1, 3, 0] = new double[] { 1, 3, 0 };
//                nodeCoords[2, 3, 0] = new double[] { 2, 3, 0 };

//                nodeCoords[0, 0, 1] = new double[] { 0, 0, 1 };
//                nodeCoords[1, 0, 1] = new double[] { 1, 0, 1 };
//                nodeCoords[2, 0, 1] = new double[] { 2, 0, 1 };
//                nodeCoords[0, 1, 1] = new double[] { 0, 1, 1 };
//                nodeCoords[1, 1, 1] = new double[] { 1, 1, 1 };
//                nodeCoords[2, 1, 1] = new double[] { 2, 1, 1 };
//                nodeCoords[0, 2, 1] = new double[] { 0, 2, 1 };
//                nodeCoords[1, 2, 1] = new double[] { 1, 2, 1 };
//                nodeCoords[2, 2, 1] = new double[] { 2, 2, 1 };
//                nodeCoords[0, 3, 1] = new double[] { 0, 3, 1 };
//                nodeCoords[1, 3, 1] = new double[] { 1, 3, 1 };
//                nodeCoords[2, 3, 1] = new double[] { 2, 3, 1 };

//                nodeCoords[0, 0, 2] = new double[] { 0, 0, 2 };
//                nodeCoords[1, 0, 2] = new double[] { 1, 0, 2 };
//                nodeCoords[2, 0, 2] = new double[] { 2, 0, 2 };
//                nodeCoords[0, 1, 2] = new double[] { 0, 1, 2 };
//                nodeCoords[1, 1, 2] = new double[] { 1, 1, 2 };
//                nodeCoords[2, 1, 2] = new double[] { 2, 1, 2 };
//                nodeCoords[0, 2, 2] = new double[] { 0, 2, 2 };
//                nodeCoords[1, 2, 2] = new double[] { 1, 2, 2 };
//                nodeCoords[2, 2, 2] = new double[] { 2, 2, 2 };
//                nodeCoords[0, 3, 2] = new double[] { 0, 3, 2 };
//                nodeCoords[1, 3, 2] = new double[] { 1, 3, 2 };
//                nodeCoords[2, 3, 2] = new double[] { 2, 3, 2 };

//                nodeCoords[0, 0, 3] = new double[] { 0, 0, 3 };
//                nodeCoords[1, 0, 3] = new double[] { 1, 0, 3 };
//                nodeCoords[2, 0, 3] = new double[] { 2, 0, 3 };
//                nodeCoords[0, 1, 3] = new double[] { 0, 1, 3 };
//                nodeCoords[1, 1, 3] = new double[] { 1, 1, 3 };
//                nodeCoords[2, 1, 3] = new double[] { 2, 1, 3 };
//                nodeCoords[0, 2, 3] = new double[] { 0, 2, 3 };
//                nodeCoords[1, 2, 3] = new double[] { 1, 2, 3 };
//                nodeCoords[2, 2, 3] = new double[] { 2, 2, 3 };
//                nodeCoords[0, 3, 3] = new double[] { 0, 3, 3 };
//                nodeCoords[1, 3, 3] = new double[] { 1, 3, 3 };
//                nodeCoords[2, 3, 3] = new double[] { 2, 3, 3 };

//                nodeCoords[0, 0, 4] = new double[] { 0, 0, 4 };
//                nodeCoords[1, 0, 4] = new double[] { 1, 0, 4 };
//                nodeCoords[2, 0, 4] = new double[] { 2, 0, 4 };
//                nodeCoords[0, 1, 4] = new double[] { 0, 1, 4 };
//                nodeCoords[1, 1, 4] = new double[] { 1, 1, 4 };
//                nodeCoords[2, 1, 4] = new double[] { 2, 1, 4 };
//                nodeCoords[0, 2, 4] = new double[] { 0, 2, 4 };
//                nodeCoords[1, 2, 4] = new double[] { 1, 2, 4 };
//                nodeCoords[2, 2, 4] = new double[] { 2, 2, 4 };
//                nodeCoords[0, 3, 4] = new double[] { 0, 3, 4 };
//                nodeCoords[1, 3, 4] = new double[] { 1, 3, 4 };
//                nodeCoords[2, 3, 4] = new double[] { 2, 3, 4 };


//                return nodeCoords[nodeIdx[0], nodeIdx[1], nodeIdx[2]];
//            }

//            public int GetNodeID(int[] nodeIdx)
//            {
//                var nodeIDs = new int[3, 4, 5];
//                nodeIDs[0, 0, 0] = 0;
//                nodeIDs[1, 0, 0] = 1;
//                nodeIDs[2, 0, 0] = 2;
//                nodeIDs[0, 1, 0] = 3;
//                nodeIDs[1, 1, 0] = 4;
//                nodeIDs[2, 1, 0] = 5;
//                nodeIDs[0, 2, 0] = 6;
//                nodeIDs[1, 2, 0] = 7;
//                nodeIDs[2, 2, 0] = 8;
//                nodeIDs[0, 3, 0] = 9;
//                nodeIDs[1, 3, 0] = 10;
//                nodeIDs[2, 3, 0] = 11;

//                nodeIDs[0, 0, 1] = 12;
//                nodeIDs[1, 0, 1] = 13;
//                nodeIDs[2, 0, 1] = 14;
//                nodeIDs[0, 1, 1] = 15;
//                nodeIDs[1, 1, 1] = 16;
//                nodeIDs[2, 1, 1] = 17;
//                nodeIDs[0, 2, 1] = 18;
//                nodeIDs[1, 2, 1] = 19;
//                nodeIDs[2, 2, 1] = 20;
//                nodeIDs[0, 3, 1] = 21;
//                nodeIDs[1, 3, 1] = 22;
//                nodeIDs[2, 3, 1] = 23;

//                nodeIDs[0, 0, 2] = 24;
//                nodeIDs[1, 0, 2] = 25;
//                nodeIDs[2, 0, 2] = 26;
//                nodeIDs[0, 1, 2] = 27;
//                nodeIDs[1, 1, 2] = 28;
//                nodeIDs[2, 1, 2] = 29;
//                nodeIDs[0, 2, 2] = 30;
//                nodeIDs[1, 2, 2] = 31;
//                nodeIDs[2, 2, 2] = 32;
//                nodeIDs[0, 3, 2] = 33;
//                nodeIDs[1, 3, 2] = 34;
//                nodeIDs[2, 3, 2] = 35;

//                nodeIDs[0, 0, 3] = 36;
//                nodeIDs[1, 0, 3] = 37;
//                nodeIDs[2, 0, 3] = 38;
//                nodeIDs[0, 1, 3] = 39;
//                nodeIDs[1, 1, 3] = 40;
//                nodeIDs[2, 1, 3] = 41;
//                nodeIDs[0, 2, 3] = 42;
//                nodeIDs[1, 2, 3] = 43;
//                nodeIDs[2, 2, 3] = 44;
//                nodeIDs[0, 3, 3] = 45;
//                nodeIDs[1, 3, 3] = 46;
//                nodeIDs[2, 3, 3] = 47;

//                nodeIDs[0, 0, 4] = 48;
//                nodeIDs[1, 0, 4] = 49;
//                nodeIDs[2, 0, 4] = 50;
//                nodeIDs[0, 1, 4] = 51;
//                nodeIDs[1, 1, 4] = 52;
//                nodeIDs[2, 1, 4] = 53;
//                nodeIDs[0, 2, 4] = 54;
//                nodeIDs[1, 2, 4] = 55;
//                nodeIDs[2, 2, 4] = 56;
//                nodeIDs[0, 3, 4] = 57;
//                nodeIDs[1, 3, 4] = 58;
//                nodeIDs[2, 3, 4] = 59;


//                return nodeIDs[nodeIdx[0], nodeIdx[1], nodeIdx[2]];
//            }

//            public int[] GetNodeIdx(int nodeID)
//            {
//                var nodeIndices = new int[,]
//                {
//                    { 0, 0, 0 },
//                    { 1, 0, 0 },
//                    { 2, 0, 0 },
//                    { 0, 1, 0 },
//                    { 1, 1, 0 },
//                    { 2, 1, 0 },
//                    { 0, 2, 0 },
//                    { 1, 2, 0 },
//                    { 2, 2, 0 },
//                    { 0, 3, 0 },
//                    { 1, 3, 0 },
//                    { 2, 3, 0 },

//                    { 0, 0, 1 },
//                    { 1, 0, 1 },
//                    { 2, 0, 1 },
//                    { 0, 1, 1 },
//                    { 1, 1, 1 },
//                    { 2, 1, 1 },
//                    { 0, 2, 1 },
//                    { 1, 2, 1 },
//                    { 2, 2, 1 },
//                    { 0, 3, 1 },
//                    { 1, 3, 1 },
//                    { 2, 3, 1 },

//                    { 0, 0, 2 },
//                    { 1, 0, 2 },
//                    { 2, 0, 2 },
//                    { 0, 1, 2 },
//                    { 1, 1, 2 },
//                    { 2, 1, 2 },
//                    { 0, 2, 2 },
//                    { 1, 2, 2 },
//                    { 2, 2, 2 },
//                    { 0, 3, 2 },
//                    { 1, 3, 2 },
//                    { 2, 3, 2 },

//                    { 0, 0, 3 },
//                    { 1, 0, 3 },
//                    { 2, 0, 3 },
//                    { 0, 1, 3 },
//                    { 1, 1, 3 },
//                    { 2, 1, 3 },
//                    { 0, 2, 3 },
//                    { 1, 2, 3 },
//                    { 2, 2, 3 },
//                    { 0, 3, 3 },
//                    { 1, 3, 3 },
//                    { 2, 3, 3 },

//                    { 0, 0, 4 },
//                    { 1, 0, 4 },
//                    { 2, 0, 4 },
//                    { 0, 1, 4 },
//                    { 1, 1, 4 },
//                    { 2, 1, 4 },
//                    { 0, 2, 4 },
//                    { 1, 2, 4 },
//                    { 2, 2, 4 },
//                    { 0, 3, 4 },
//                    { 1, 3, 4 },
//                    { 2, 3, 4 }
//                };

//                return new int[] { nodeIndices[nodeID, 0], nodeIndices[nodeID, 1], nodeIndices[nodeID, 2] };
//            }
//        }
//    }
//}
