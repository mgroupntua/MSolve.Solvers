//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Net;
//using System.Text;
//using MGroup.MSolve.Discretization.Mesh;
//using MGroup.Solvers.DDM.Mesh;
//using MGroup.Solvers.DDM.Output;
//using MGroup.Solvers.DDM.Output.VTK;
//using Xunit;

//namespace MGroup.Solvers.DDM.Tests.Mesh
//{
//    public static class UniformMesh2DTests
//    {
//        //[Fact]
//        public static void PlotMesh()
//        {
//            var writer = new VtkMeshWriter();

//            string path = @"C:\Users\Serafeim\Desktop\PFETIDP\meshes\mesh2D.vtk";
//            double[] minCoords = { 0, 0 };
//            double[] maxCoords = { 12, 12 };
//            int[] numElements = { 3, 4 };
//            var mesh = new UniformMesh2D.Builder(minCoords, maxCoords, numElements).SetMajorAxis(1).BuildMesh();
//            //var mesh = new UniformMesh2D.Builder(minCoords, maxCoords, numElements).BuildMesh();
//            writer.WriteMesh(path, mesh, 2);
//        }

//        [Fact]
//        public static void TestEnumerateElements()
//        {
//            var mockMesh = new MockMesh2x3();
//            var mesh = new UniformMesh2D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
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
//            var mockMesh = new MockMesh2x3();
//            var mesh = new UniformMesh2D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
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
//                for (int d = 0; d < 2; ++d)
//                {
//                    Assert.Equal(coordsExpected[d], coords[d]);
//                }
//            }
//        }

//        [Fact]
//        public static void TestGetNodeID()
//        {
//            var mockMesh = new MockMesh2x3();
//            var mesh = new UniformMesh2D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
//                .BuildMesh();

//            for (int i = 0; i < mesh.NumNodes[0]; ++i)
//            {
//                for (int j = 0; j < mesh.NumNodes[1]; ++j)
//                {
//                    var idx = new int[] { i, j };
//                    int expected = mockMesh.GetNodeID(idx);
//                    int computed = mesh.GetNodeID(idx);
//                    Assert.Equal(expected, computed);
//                }
//            }
//        }

//        [Fact]
//        public static void TestGetNodeIdx()
//        {
//            var mockMesh = new MockMesh2x3();
//            var mesh = new UniformMesh2D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
//                .BuildMesh();

//            for (int id = 0; id < mesh.NumNodes[0] * mesh.NumNodes[1]; ++id)
//            {
//                int[] expected = mockMesh.GetNodeIdx(id);
//                int[] computed = mesh.GetNodeIdx(id);
//                for (int d = 0; d < 2; ++d)
//                {
//                    Assert.Equal(expected[d], computed[d]);
//                }
//            }
//        }

//        [Fact]
//        public static void TestGetNodeCoordinates()
//        {
//            var mockMesh = new MockMesh2x3();
//            var mesh = new UniformMesh2D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
//                .BuildMesh();

//            for (int i = 0; i < mesh.NumNodes[0]; ++i)
//            {
//                for (int j = 0; j < mesh.NumNodes[1]; ++j)
//                {
//                    var idx = new int[] { i, j };
//                    double[] expected = mockMesh.GetNodeCoordinates(idx);
//                    double[] computed = mesh.GetNodeCoordinates(idx);
//                    for (int d = 0; d < 2; ++d)
//                    {
//                        Assert.Equal(expected[d], computed[d]);
//                    }
//                }
//            }
//        }

//        [Fact]
//        public static void TestGetElementID()
//        {
//            var mockMesh = new MockMesh2x3();
//            var mesh = new UniformMesh2D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
//                .BuildMesh();

//            for (int i = 0; i < mesh.NumElements[0]; ++i)
//            {
//                for (int j = 0; j < mesh.NumElements[1]; ++j)
//                {
//                    var idx = new int[] { i, j };
//                    int expected = mockMesh.GetElementID(idx);
//                    int computed = mesh.GetElementID(idx);
//                    Assert.Equal(expected, computed);
//                }
//            }
//        }

//        [Fact]
//        public static void TestGetElementIdx()
//        {
//            var mockMesh = new MockMesh2x3();
//            var mesh = new UniformMesh2D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
//                .BuildMesh();

//            for (int id = 0; id < mesh.NumElements[0] * mesh.NumElements[1]; ++id)
//            {
//                int[] expected = mockMesh.GetElementIdx(id);
//                int[] computed = mesh.GetElementIdx(id);
//                for (int d = 0; d < 2; ++d)
//                {
//                    Assert.Equal(expected[d], computed[d]);
//                }
//            }
//        }

//        [Fact]
//        public static void TestGetElementConductivity()
//        {
//            var mockMesh = new MockMesh2x3();
//            var mesh = new UniformMesh2D.Builder(mockMesh.MinCoordinates, mockMesh.MaxCoordinates, mockMesh.NumElements)
//                .BuildMesh();

//            for (int i = 0; i < mesh.NumElements[0]; ++i)
//            {
//                for (int j = 0; j < mesh.NumElements[1]; ++j)
//                {
//                    var idx = new int[] { i, j };
//                    int[] expected = mockMesh.GetElementConnectivity(idx);
//                    int[] computed = mesh.GetElementConnectivity(idx);
//                    for (int d = 0; d < expected.Length; ++d)
//                    {
//                        Assert.Equal(expected[d], computed[d]);
//                    }
//                }
//            }
//        }

//        private class MockMesh2x3 : IStructuredMesh
//        {
//            public CellType CellType => CellType.Quad4;

//            public double[] MinCoordinates => new double[] { 0.0, 0.0 };

//            public double[] MaxCoordinates => new double[] { 2.0, 3.0 };

//            public int[] NumElements => new int[] { 2, 3 };

//            public int[] NumNodes => new int[] { 3, 4 };

//            public int NumElementsTotal => 6;

//            public int NumNodesTotal => 12;

//            public int NumNodesPerElement => 4;

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
//                var elementNodes = new int[,][]
//                {
//                    { new int[] { 0, 1, 4, 3 }, new int[] { 3, 4, 7, 6 }, new int[] { 6, 7, 10, 9 } },
//                    { new int[] { 1, 2, 5, 4 }, new int[] { 4, 5, 8, 7 }, new int[] { 7, 8, 11, 10 } }
//                };

//                return elementNodes[elementIdx[0], elementIdx[1]];
//            }

//            public int GetElementID(int[] elementIdx)
//            {
//                var elemIDs = new int[,]
//                {
//                    { 0, 2, 4 },
//                    { 1, 3, 5 }
//                };

//                return elemIDs[elementIdx[0], elementIdx[1]];
//            }

//            public int[] GetElementIdx(int elementID)
//            {
//                var elemIndices = new int[,]
//                {
//                    { 0, 0 },
//                    { 1, 0 },
//                    { 0, 1 },
//                    { 1, 1 },
//                    { 0, 2 },
//                    { 1, 2 },
//                };

//                return new int[] { elemIndices[elementID, 0], elemIndices[elementID, 1] };
//            }

//            public double[] GetNodeCoordinates(int[] nodeIdx)
//            {
//                var nodeCoords = new double[,][]
//                {
//                    { new double[] { 0, 0 }, new double[] { 0, 1 }, new double[] { 0, 2 }, new double[] { 0, 3 } },
//                    { new double[] { 1, 0 }, new double[] { 1, 1 }, new double[] { 1, 2 }, new double[] { 1, 3 } },
//                    { new double[] { 2, 0 }, new double[] { 2, 1 }, new double[] { 2, 2 }, new double[] { 2, 3 } },
//                };

//                return nodeCoords[nodeIdx[0], nodeIdx[1]];
//            }

//            public int GetNodeID(int[] nodeIdx)
//            {
//                var nodeIDs = new int[,]
//                {
//                    { 0, 3, 6, 9 },
//                    { 1, 4, 7, 10 },
//                    { 2, 5, 8, 11 }
//                };

//                return nodeIDs[nodeIdx[0], nodeIdx[1]];
//            }

//            public int[] GetNodeIdx(int nodeID)
//            {
//                var nodeIndices = new int[,]
//                {
//                    { 0, 0 },
//                    { 1, 0 },
//                    { 2, 0 },
//                    { 0, 1 },
//                    { 1, 1 },
//                    { 2, 1 },
//                    { 0, 2 },
//                    { 1, 2 },
//                    { 2, 2 },
//                    { 0, 3 },
//                    { 1, 3 },
//                    { 2, 3 },
//                };

//                return new int[] { nodeIndices[nodeID, 0], nodeIndices[nodeID, 1] };
//            }
//        }
//    }
//}
