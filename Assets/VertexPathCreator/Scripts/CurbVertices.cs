using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;
using PathCreation;

public class CurbVertices : MonoBehaviour
{
    public PathCreator pathCreator;

    protected VertexPath path => pathCreator.path;

    /*
 * Mesh To Extrude
 * Mesh which will be extruded along the spline
 */
    [SerializeField]
    private Mesh _meshToExtrude;
    public Mesh MeshToExtrude
    { get { return _meshToExtrude; } set { _meshToExtrude = value; } }

    /*
     * Curb Mesh
     */
    [SerializeField]
    private Mesh curbMesh;
    public Mesh CurbMesh
    { get { return curbMesh; } set { curbMesh = value; } }

    /*
    * Material which is applied to the curb mesh 
    */
    [SerializeField]
    private Material _curbMaterial;
    public Material CurbMaterial
    { get { return _curbMaterial; } set { _curbMaterial = value; } }

    /*
     * Store a GameObject for the track
     */
    [SerializeField]
    private Mesh curbsMesh;

    /*
     * Class to store all the procedural data
     * Vertices
     * Triangles
     * UVs
     */
    ProMesh proMesh;

    [SerializeField]
    private GameObject curb;

    /// <summary>
    /// Create the track
    /// </summary>
    /// 

    [SerializeField]
    bool CreateTrackOnStar;
    [SerializeField]
    bool onPlanXZ;

    private void Start()
    {
        if (CreateTrackOnStar)
            CreateTrack();
    }

    [ContextMenu("Create Mesh Track")]
    public void CreateTrack()
    {
        if (path.NumPoints > 0)
        {
            proMesh = new ProMesh();
            ProTwoDShape pTwoDShape = CreateTwoDShape(MeshToExtrude);
            CreateTrackVertices();
            CreateTackMesh(pTwoDShape);
            CreateCurbVertices();
        }
        else
        {
            Debug.Log("SPLICE CREATOR: There is no spline. Make sure there is a Voronoi diagram gnerated. Then press the" +
                " 'Generate Spline'");
        }
    }

    /// <summary>
    /// Create a new ProTwoDShape Object to store 
    /// all the data for the mesh to extrude along the spline 
    /// </summary>
    /// <returns></returns>
    private ProTwoDShape CreateTwoDShape(Mesh a_mesh)
    {
        //ref to futhest vertex on the left and right
        Vector2 farLeftV = Vector2.zero;
        Vector2 farRightV = Vector2.zero;
        //Create List of all vertices
        List<Vector2> twoDShape = new List<Vector2>();

        for (int i = 0; i < a_mesh.vertexCount;)
        {
            //Add each vertex to twoDShape List
            twoDShape.Add(new Vector2(a_mesh.vertices[i].x,
                                      a_mesh.vertices[i].y));
            //Check if current vertex is further left then any other
            if (twoDShape[twoDShape.Count - 1].x < farLeftV.x)
            {
                farLeftV = twoDShape[twoDShape.Count - 1];
            }
            //Check if current vertex is further right then any other
            if (twoDShape[twoDShape.Count - 1].x > farRightV.x)
            {
                farRightV = twoDShape[twoDShape.Count - 1];
            }

            //if we are on the first vertex
            //then inceremnt i by 3
            if (i == 0)
            {
                i += 3;
            }
            //otherwist incerment i by 2
            //this done due to how the vertices are layed out
            //
            //1     2       4
            //
            //0     3       5
            else
            {
                i += 2;
            }
        }

        //Get the overall width of the 2d shape
        float twoDShapeWidth = Mathf.Abs(farLeftV.x - farRightV.x);

        return new ProTwoDShape(twoDShape, twoDShapeWidth, (farLeftV.x - farRightV.x));
    }

    /// <summary>
    /// Create all the vertices for the mesh
    /// </summary>
    private void CreateTrackVertices()
    {
        ProTwoDShape pTwoDShape = CreateTwoDShape(MeshToExtrude);

        //Create an array of OrientedPoint
        OrientedPoint[] oPoints = new OrientedPoint[path.NumPoints];
        //Assign each oPoint 
        for (int i = 0; i < oPoints.Length; i++)
        {
            //Get the point before current point
            Vector3 lastPoint = path.GetPoint(((i - 1) + oPoints.Length) % oPoints.Length);
            //Vector3 lastPoint = _spline.AllSmoothPoints[];
            //Get current point
            Vector3 currentPoint = path.GetPoint(i);
            //Vector3 currentPoint = _spline.AllSmoothPoints[i];
            //Get the next point 
            Vector3 nextPoint = path.GetPoint(((i + 1) + oPoints.Length) % oPoints.Length);
            //Vector3 nextPoint = _spline.AllSmoothPoints[];

            //Get the direction from currentPoint to lastPoint
            Vector3 lc = currentPoint - lastPoint;
            //Get the direction from nextPoint to currentPoint
            Vector3 cn = nextPoint - currentPoint;

            //Avgrage the start and end points for each curve
            //Set the points direction. This direction is the forward vector 
            Vector3 dir = (lc + cn) * 0.5f;
            dir.Normalize();

            //Get the right Vector 
            dir = Vector3.Cross(dir, Vector3.up).normalized;

            //new the oPoint
            oPoints[i] = new OrientedPoint(path.GetPoint(i),
                                           Quaternion.LookRotation(dir, Vector3.up));
        }

        //Generate all the center points on the mesh
        //This will also set the forward rotation
        //and UV for each point
        for (int i = 0; i < oPoints.Length; i++)
        {
            for (int j = 0; j < pTwoDShape.Vertices.Count; j++)
            {
                //get the forward vector
                Vector3 fwd = oPoints[i].Rotation * Vector3.forward;
                //get the position
                Vector3 pos = oPoints[i].Position + (pTwoDShape.Vertices[j].x * fwd);
                pos.y = oPoints[i].Position.y + pTwoDShape.Vertices[j].y;

                //add the new vertex to the proMesh.Vertices list
                proMesh.Vertices.Add(pos);

                //Setup the uv coords for this vertex
                float x = (float)j / (pTwoDShape.Vertices.Count - 1);
                float y = (float)i / oPoints.Length;
                Vector2 uv = new Vector2(x, y);
                proMesh.Uvs.Add(uv);
            }
        }
    }

    /// <summary>
    /// Create the mesh 
    /// </summary>
    /// <param name="a_vertices"></param>
    /// <param name="a_shapeToExtrude"></param>
    /// <returns></returns>
    private void CreateTackMesh(ProTwoDShape a_twoDShape)
    {
        //loop though all the vertices
        for (int i = 0; i < proMesh.Vertices.Count; i++)
        {
            //this is to check if we are on the last vertex of the twoDShape
            //e.x ((1 + 1) / 2) % 1 = 1
            //we are on the end of the twoDShape
            //we can not add trianles
            //so continue
            float value = ((i + 1) / (float)a_twoDShape.Vertices.Count) % 1f;
            if (i != 0 && value == 0)
            {
                continue;
            }

            //Make tri
            proMesh.Triangles.Add(i);
            proMesh.Triangles.Add(i + a_twoDShape.Vertices.Count);
            proMesh.Triangles.Add(i + a_twoDShape.Vertices.Count + 1);

            proMesh.Triangles.Add(i);
            proMesh.Triangles.Add(i + a_twoDShape.Vertices.Count + 1);
            proMesh.Triangles.Add(i + 1);
        }
    }


    /// <summary>
    /// Create all the vertices for the curbs around the track
    /// </summary>
    private void CreateCurbVertices()
    {
        ProTwoDShape curbMesh = CreateTwoDShape(CurbMesh);

        float center = curbMesh.Center;

        for (int i = 0; i < curbMesh.Vertices.Count; i++)
        {
            Vector2 v = curbMesh.Vertices[i];
            v *= 100;
            curbMesh.Vertices[i] = v;
        }
        curbMesh.ShapeWidth *= 100;

        //Store all the curb meshes made
        List<CombineInstance> allCurbMeshes = new List<CombineInstance>();
        if (true)
        {
            //Go though all the vertices in the track mesh
            if (true)
            {
                int i = 0;
                //Store all the vertices for the current curb
                List<Vector3> curbVertices = new List<Vector3>();
                //Store all the uvs for the current curb
                List<Vector2> curbUvs = new List<Vector2>();
                //Store all the triangle for the current curb
                List<int> curbTriangles = new List<int>();

                //Get the number of segmnets each spline curve it broken
                //up into
                int segStart = proMesh.Vertices.Count;

                //if we are at the start of a curve 
                if (i % segStart == 0)
                {
                    //Get the distance for the left side of the spline
                    float splineLeft = (proMesh.Vertices[i] -
                                        proMesh.Vertices[Helper.LoopIndex(i, segStart - 2, proMesh.Vertices.Count)]).magnitude;
                    //Get the distance for the right side of the spline
                    float splineRight = (proMesh.Vertices[Helper.ALoopIndex(i, proMesh.Vertices.Count)] -
                                        proMesh.Vertices[Helper.LoopIndex(i, segStart - 1, proMesh.Vertices.Count)]).magnitude;

                    //check which side of the track mesh is smaller
                    //this is the side which will have the curb made for
                    if (true)
                    {
                        curbMesh.Vertices.Reverse();
                        //Setup the vertices
                        for (int j = 0; j < segStart ; j += 2)
                        {
                            for (int k = 0; k < curbMesh.Vertices.Count; k++)
                            {
                                Vector3 cV;
                                Vector3 cVT;
                                Vector3 dir;
                                if (j == 0 || j == segStart - 2)
                                {
                                     cV = proMesh.Vertices[Helper.LoopIndex(i, j + 1, proMesh.Vertices.Count)];
                                     cVT = proMesh.Vertices[Helper.LoopIndex(i, j, proMesh.Vertices.Count)];
                                }
                                else
                                {
                                    cV = proMesh.Vertices[Helper.LoopIndex(i, j, proMesh.Vertices.Count)];
                                    cVT = proMesh.Vertices[Helper.LoopIndex(i, j + 1, proMesh.Vertices.Count)];
                                }

                                dir = cV - cVT;
                                //ROTATE VERTICES

                                //Set the new vertex position
                                Vector3 pos;
                                if (curbMesh.Vertices[k].x < 0)
                                    pos = cVT + (curbMesh.Vertices[k].x * dir.normalized);
                                else
                                    pos = cV + (curbMesh.Vertices[k].x * dir.normalized);

                                if(onPlanXZ)
                                    pos.y = 0 + curbMesh.Vertices[k].y;
                                else
                                    pos.y += curbMesh.Vertices[k].y;

                                //add curb to vertex list
                                curbVertices.Add(pos);

                                //Set the new vertex uv
                                float x = (float)j / (curbMesh.Vertices.Count - 1);
                                float y = (float)i / 4 + 1;
                                Vector2 uv = new Vector2(x, y);
                                curbUvs.Add(uv);
                            }
                        }
                        //for (int j = segStart - 2; j < segStart; j += 2)
                        //{
                        //    for (int k = 0; k < curbMesh.Vertices.Count; k++)
                        //    {
                        //        //ROTATE VERTICES
                        //        Vector3 cV = proMesh.Vertices[Helper.LoopIndex(i, j + 1, proMesh.Vertices.Count)];
                        //        Vector3 cVT = proMesh.Vertices[Helper.LoopIndex(i, j, proMesh.Vertices.Count)];
                        //        Vector3 dir = cV -
                        //                      cVT;

                        //        //Set the new vertex position
                        //        Vector3 pos;
                        //        if (curbMesh.Vertices[k].x < 0)
                        //            pos = cVT + (curbMesh.Vertices[k].x * dir.normalized);
                        //        else
                        //            pos = cV + (curbMesh.Vertices[k].x * dir.normalized);


                        //        pos.y = 0 + curbMesh.Vertices[k].y;


                        //        //add curb to vertex list
                        //        curbVertices.Add(pos);

                        //        //Set the new vertex uv
                        //        float x = (float)j / (curbMesh.Vertices.Count - 1);
                        //        float y = (float)i / 4 + 1;
                        //        Vector2 uv = new Vector2(x, y);
                        //        curbUvs.Add(uv);
                        //    }
                        //}
                        //Setup the triangles
                        for (int j = 0; j < curbVertices.Count - (curbMesh.Vertices.Count); j++)
                        {
                            float value = (j + 1) / (float)curbMesh.Vertices.Count % 1f;
                            if (value != 0)
                            {
                                curbTriangles.Add(j);
                                curbTriangles.Add(j + curbMesh.Vertices.Count);
                                curbTriangles.Add(j + curbMesh.Vertices.Count + 1);
                                curbTriangles.Add(j);
                                curbTriangles.Add(j + curbMesh.Vertices.Count + 1);
                                curbTriangles.Add(j + 1);
                            }
                        }
                        curbMesh.Vertices.Reverse();
                    }

                    CombineInstance ci = new CombineInstance();
                    Mesh m = new Mesh();
                    m.vertices = curbVertices.ToArray();
                    m.uv = curbUvs.ToArray();
                    m.triangles = curbTriangles.ToArray();
                    ci.mesh = m;
                    ci.transform = transform.localToWorldMatrix;
                    allCurbMeshes.Add(ci);
                }
            }
        }

        curbsMesh = new Mesh();
        //Combine all curb meshs in to one
        curbsMesh.CombineMeshes(allCurbMeshes.ToArray());

        //Clean up all the meshes by destroying them
        for (int i = 0; i < allCurbMeshes.Count; i++)
        {
            allCurbMeshes[i].mesh.Clear();
        }

        //then create a new GameObject add a mesh filter and mesh renderer to it
        if (!curb)
            curb = new GameObject("Curbs");
        if (!curb.GetComponent<MeshFilter>())
            curb.AddComponent<MeshFilter>();
        if (!curb.GetComponent<MeshRenderer>())
            curb.AddComponent<MeshRenderer>();

        curbsMesh.RecalculateNormals();

        //set the mesh filter and mesh renderer
        curb.GetComponent<MeshFilter>().sharedMesh = curbsMesh;
        curb.GetComponent<MeshRenderer>().sharedMaterial = CurbMaterial;
    }

    /// <summary>
    /// Check if any there are any overlapping
    /// </summary>
    private void CheckForOverLap()
    {
        //set splineIntersects
        bool splineIntersects = false;
        //set splineIntersectsStart
        int splineIntersectsStart = -1;
        //create a new list of ints
        List<int> pointsToChange = new List<int>();

        for (int i = 3; i < proMesh.Vertices.Count; i += 2)
        {
            int nLineIndex = ((i + 1) + proMesh.Vertices.Count) % proMesh.Vertices.Count;
            int nLineIndexPlusOne = ((i + 2) + proMesh.Vertices.Count) % proMesh.Vertices.Count;

            //get the before, current and next line
            Line bLine = new Line
            {
                p0 = new Vector2(proMesh.Vertices[i - 3].x, proMesh.Vertices[i - 3].z),
                p1 = new Vector2(proMesh.Vertices[i - 2].x, proMesh.Vertices[i - 2].z)
            };
            Line cLine = new Line
            {
                p0 = new Vector2(proMesh.Vertices[i - 1].x, proMesh.Vertices[i - 1].z),
                p1 = new Vector2(proMesh.Vertices[i].x, proMesh.Vertices[i].z)
            };
            Line nLine = new Line
            {
                p0 = new Vector2(proMesh.Vertices[nLineIndex].x, proMesh.Vertices[nLineIndex].z),
                p1 = new Vector2(proMesh.Vertices[nLineIndexPlusOne].x, proMesh.Vertices[nLineIndexPlusOne].z)
            };

            Vector2 intersects;
            //check for intersects from before and current line
            if (Helper.LineInsterection(bLine.p0, bLine.p1, cLine.p0, cLine.p1, out intersects))
            {
                //set splineIntersects to true
                splineIntersects = true;
                //check if splineIntersectsStart is -1
                //set splineIntersectsStart
                if (splineIntersectsStart == -1)
                {
                    splineIntersectsStart = i;
                }
                //add two points to pointsToChange
                pointsToChange.Add(i - 1);
                pointsToChange.Add(i - 3);

                //set newPos
                Vector3 newPos = (proMesh.Vertices[((i - 5) + proMesh.Vertices.Count) % proMesh.Vertices.Count] +
                                 proMesh.Vertices[nLineIndex]) * 0.5f;
            }
            else
            {
                //if splineIntersects is true 
                if (splineIntersects)
                {
                    //i - splineIntersectsStart is greater than 6 then we are on a 
                    //ne curve
                    if (i - splineIntersectsStart > 6)
                    {
                        //set splineIntersects
                        splineIntersects = false;
                        //set splineIntersectsStart
                        splineIntersectsStart = -1;
                        //create a vector3 for the center point
                        Vector3 center = Vector3.zero;

                        //add all the pointsToChange to center
                        for (int k = 0; k < pointsToChange.Count; k++)
                        {
                            center += proMesh.Vertices[pointsToChange[k]];
                        }
                        //divide center by the number of pointsToChange
                        center /= pointsToChange.Count;
                        //set all the pointsToChange to the center
                        for (int k = 0; k < pointsToChange.Count; k++)
                        {
                            proMesh.Vertices[pointsToChange[k]] = center;
                        }
                    }
                    //clear pointsToChange
                    pointsToChange.Clear();
                }
            }
        }
    }

    /// <summary>
    /// return a new point in a bezier curve with 3 points
    /// </summary>
    /// <param name="a_P0"></param>
    /// <param name="a_P1"></param>
    /// <param name="a_P2"></param>
    /// <param name="a_t"></param>
    /// <returns></returns>
    private Vector3 GetPointInSmoothSpline(Vector3 a_P0, Vector3 a_P1, Vector3 a_P2, float a_t)
    {
        //(1 - a_t)2 P0 + 2(1-t)tP1 + tP2
        //float u = 1 - a_t;
        //float tt = a_t * a_t;
        //float uu = u - u;
        //Vector3 p = uu * a_P0;
        //p += 2 * u * a_t * a_P1;
        //p += tt * a_P2;
        //return p;
        return Vector3.Lerp(Vector3.Lerp(a_P0, a_P1, a_t), Vector3.Lerp(a_P1, a_P2, a_t), a_t);
    }
}

/// <summary>
/// Class for the 2d shape
/// This holds all the vertices and the width of the shape
/// </summary>
public class ProTwoDShape
{
    public List<Vector2> Vertices;
    public float ShapeWidth;
    public float Center;

    /// <summary>
    /// Constructor
    /// </summary>
    /// <param name="a_vertices"></param>
    /// <param name="a_width"></param>
    public ProTwoDShape(List<Vector2> a_vertices, float a_width, float center = 0)
    {
        Vertices = a_vertices;
        ShapeWidth = a_width;
        Center = center;
    }
}

/// <summary>
/// Pro Mesh Object. This object holds all data need for a Mesh
/// -Vertices
/// -Triangles
/// -Uvs
/// </summary>
public class ProMesh
{
    public List<Vector3> Vertices;
    public List<int> Triangles;
    public List<Vector2> Uvs;

    public ProMesh()
    {
        Vertices = new List<Vector3>();
        Triangles = new List<int>();
        Uvs = new List<Vector2>();
    }
}

/// <summary>
/// OrientedPoint class. This class stores a position and a rotation
/// </summary>
public class OrientedPoint
{
    public Vector3 Position;
    public Quaternion Rotation;

    public OrientedPoint(Vector3 a_position, Quaternion a_rotation)
    {
        Position = a_position;
        Rotation = a_rotation;
    }

    public OrientedPoint(Vector3 a_position, Vector3 a_rotation)
    {
        Position = a_position;
        Rotation = Quaternion.Euler(a_rotation.x, a_rotation.y, a_rotation.z);
    }
}

/// <summary>
/// Struct for a line object
/// This contains two points. Both are Vector2s
/// </summary>
public struct Line
{
    public Vector2 p0;
    public Vector2 p1;

    public Vector3 P0()
    {
        return new Vector3(p0.x, 0, p0.y);
    }
    public Vector3 P1()
    {
        return new Vector3(p1.x, 0, p1.y);
    }
}