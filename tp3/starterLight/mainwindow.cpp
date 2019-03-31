#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>

float max_gauss = 0;
float min_gauss = 0;
float moy = 0;

/* **** début de la partie à compléter **** */
float MainWindow::faceArea(MyMesh* _mesh, int faceID)
{
    QVector<float> vectors;
    float area = 0.0;

    FaceHandle fh = _mesh->face_handle(faceID);

    for(MyMesh::FaceVertexIter curVertex = _mesh->fv_iter(fh); curVertex.is_valid(); curVertex ++)
    {
        VertexHandle vh = *curVertex;
        vectors.append(_mesh->point(vh)[0]);
        vectors.append(_mesh->point(vh)[1]);
        vectors.append(_mesh->point(vh)[2]);
    }
    //on a nos 3 sommets de chaque face : vertex a = vector[0][1][2] vertex b = [3][4][5] vertex c = [6][7][8]

    QVector<float> vectAB;
    QVector<float> vectAC;
    QVector<float> produitVect;

    //On calcul les vecteur AB et AC
    for(int i=0; i<3; i++)
    {
        vectAB.append(vectors[i+3]-vectors[i]);
        vectAC.append(vectors[i+6]-vectors[i]);
    }
    // AB : (x1,y1,z1)    AC : (x2,y2,z2)
    produitVect.append(vectAB[1]*vectAC[2]-vectAB[2]*vectAC[1]); //produitVectoriel : vx
    produitVect.append(vectAB[0]*vectAC[2]-vectAB[2]*vectAC[0]); //produitVectoriel: vy
    produitVect.append(vectAB[0]*vectAC[1]-vectAB[1]*vectAC[0]); //produitVectoriel: vz

    //calcul area = 1/2||produiVect||
    area = 0.5*sqrt(produitVect[0]*produitVect[0]+produitVect[1]*produitVect[1]+produitVect[2]*produitVect[2]);

    return area;
}

float MainWindow::barycentriqueArea(MyMesh* _mesh, int vertexID)
{
    float baryarea = 0.0;
    float tierce = (float)1/(float)3;
    VertexHandle vert = _mesh->vertex_handle(vertexID);
    //on itere sur toute les face du sommet envoyé en parametre
    for(MyMesh::VertexFaceIter vf = _mesh->vf_iter(vert); vf.is_valid(); vf++)
    {
        FaceHandle currentFace = *vf;
        baryarea += tierce*(faceArea(_mesh,currentFace.idx()));
    }
    return baryarea;
}

float MainWindow::angleFF(MyMesh* _mesh, int faceID0,  int faceID1, int vertID0, int vertID1)
{
    /* **** à compléter ! **** */
    return 0.0;
}

float MainWindow::angleEE(MyMesh* _mesh, int vertexID,  int faceID)
{

    VertexHandle vh = _mesh->vertex_handle(vertexID);
    FaceHandle fh = _mesh->face_handle(faceID);
    QVector<VertexHandle> listePoints;
    QVector<VertexHandle> listePointsOnFace;
    QVector<float> vectors;
    //On identifie les point voisins de vertexID qui appartiennent à faceID

    //tout les points voisins de vertexID
    for (MyMesh::VertexVertexIter curVertex = _mesh->vv_iter(vh); curVertex.is_valid(); curVertex ++)
    {
        VertexHandle v = *curVertex;
        listePoints.append(v);
    }
    //parmis ces points ceux qui appartiennent à la faceID

    for(MyMesh::FaceVertexIter curVertex = _mesh->fv_iter(fh); curVertex.is_valid(); curVertex ++)
    {
        VertexHandle v = *curVertex;
        if(listePoints.contains(v) && v.idx() != vertexID)
        {
            listePointsOnFace.append(v);
        }
    }

    //On créer des vecteurs a partir des point obtenu
    for(int i=0; i<listePointsOnFace.size();i++)
    {
        vectors.append((_mesh->point(listePointsOnFace[i])[0])-(_mesh->point(vh)[0]));
        vectors.append((_mesh->point(listePointsOnFace[i])[1])-(_mesh->point(vh)[1]));
        vectors.append((_mesh->point(listePointsOnFace[i])[2])-(_mesh->point(vh)[2]));

    }

    //on normalise les vecteurs obtenu
    QVector<float> norme;

    float tmp = vectors[0]*vectors[0]+vectors[1]*vectors[1]+vectors[2]*vectors[2];
    norme.append(sqrt(tmp)); //sqrt(x^2+y^2+z^2) : C'est la norme du premier vecteur

    tmp = vectors[3]*vectors[3]+vectors[4]*vectors[4]+vectors[5]*vectors[5];
    norme.append(sqrt(tmp)); //C'est la norme du second vecteurs

    //normalisation des vecteur :  chacun des vecteur / par la norme
    for(int i=0; i<listePointsOnFace.size(); i++) //2 vecteur : 2 iterations : 6points
    {
        vectors[i*3] = vectors[i*3]/norme[i];
        vectors[i*3+1] = vectors[i*3+1]/norme[i];
        vectors[i*3+2] = vectors[i*3+2]/norme[i];
    }

    float a = vectors[0]*vectors[0+3];
    float b = vectors[1]*vectors[1+3];
    float c = vectors[2]*vectors[2+3];

    float prodScal = a+b+c;
    abs(prodScal); //garantie que le produit scalaire soit positif

    float angle = acos(prodScal)/**180/M_PI*/;
    return angle; //en radians
}

void MainWindow::H_Curv(MyMesh* _mesh)
{
    //Non completé
}

void MainWindow::K_Curv(MyMesh* _mesh)
{
    float interval = val_med(_mesh);
    resetAllColorsAndThickness(_mesh);

    float sommeAngle = 0;
    float gaussFace_Value;

    //on doit iterer sur chacun des sommet de toutes les face

    for(MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        FaceHandle fh = *curFace;
        for (MyMesh::FaceVertexIter curVertex = _mesh->fv_iter(fh); curVertex.is_valid(); curVertex ++)
        {
            VertexHandle vh = *curVertex;
            sommeAngle += angleEE(_mesh, vh.idx(), fh.idx());

        }

        for(MyMesh::FaceVertexIter curVertex = _mesh->fv_iter(fh); curVertex.is_valid(); curVertex ++)
        {
            VertexHandle vh2 = *curVertex;

            float gaussValue = (float(1)/barycentriqueArea(_mesh,vh2.idx()))*(2*M_PI-sommeAngle);

            gaussFace_Value += gaussValue;

        }

        gaussFace_Value = (gaussFace_Value/3)*255/((max_gauss)); // Produit en croix
        qDebug()<<gaussFace_Value<<endl;

        if (gaussFace_Value <0)
             gaussFace_Value =0;
        if (gaussFace_Value > 255)
            gaussFace_Value = 255;

         _mesh->set_color(fh, MyMesh::Color(gaussFace_Value/2, gaussFace_Value/4, gaussFace_Value/3));

        sommeAngle = 0;
    }

    //qDebug() << max_gauss <<endl;

    displayMesh(_mesh);
}

float MainWindow::val_med(MyMesh *_mesh)
{
    float sommeAngle = 0;
    float gaussFace_Value;
    for(MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        FaceHandle fh = *curFace;
        for (MyMesh::FaceVertexIter curVertex = _mesh->fv_iter(fh); curVertex.is_valid(); curVertex ++)
        {
            VertexHandle vh = *curVertex;
            sommeAngle += angleEE(_mesh, vh.idx(), fh.idx());

        }

        for(MyMesh::FaceVertexIter curVertex = _mesh->fv_iter(fh); curVertex.is_valid(); curVertex ++)
        {
            VertexHandle vh2 = *curVertex;

            float gaussValue = (float(1)/barycentriqueArea(_mesh,vh2.idx()))*(2*M_PI-sommeAngle);
            //GaussArray.append(gaussValue);


            if(max_gauss < gaussValue)
            {
                 max_gauss = gaussValue;
            }
            if(min_gauss > gaussValue)
            {
                min_gauss = gaussValue;
            }

            gaussFace_Value += gaussValue;



        }


      }

    moy = abs(static_cast <float>(gaussFace_Value/_mesh->n_faces()));
    return abs(static_cast <float>(max_gauss/_mesh->n_faces()));


}
/* **** fin de la partie à compléter **** */



/* **** début de la partie boutons et IHM **** */
void MainWindow::on_pushButton_H_clicked()
{
    //test perso
    //qDebug()<<faceArea(&mesh,0)<<endl;
    //qDebug()<<barycentriqueArea(&mesh,0)<<endl;
    //qDebug()<<"angle"<<angleEE(&mesh,0,0)<<endl;


    //fin de test

    H_Curv(&mesh);
    displayMesh(&mesh, true); // true permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}

void MainWindow::on_pushButton_K_clicked()
{
    K_Curv(&mesh);
    //displayMesh(&mesh, true); // true permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}

/*
    Cette fonction est à utiliser UNIQUEMENT avec le fichier testAngleArea.obj
    Elle est appelée par le bouton "Test angles/aires"

    Elle permet de vérifier les fonctions faceArea, angleFF et angleEE.
    Elle doit afficher :

    Aire de la face 0 : 2
    Aire de la face 1 : 2
    Angle entre les faces 0 et 1 : 1.5708
    Angle entre les faces 1 et 0 : -1.5708
    Angle au sommet 1 sur la face 0 : 0.785398
*/

void MainWindow::on_pushButton_angleArea_clicked()
{
    qDebug() << "Aire de la face 0 :" << faceArea(&mesh, 0);
    qDebug() << "Aire de la face 1 :" << faceArea(&mesh, 1);

    qDebug() << "Angle entre les faces 0 et 1 :" << angleFF(&mesh, 0, 1, 1, 2);
    qDebug() << "Angle entre les faces 1 et 0 :" << angleFF(&mesh, 1, 0, 1, 2);

    qDebug() << "Angle au sommet 1 sur la face 0 :" << angleEE(&mesh, 1, 0);
    qDebug() << "Angle au sommet 3 sur la face 1 :" << angleEE(&mesh, 3, 1);
}

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}
/* **** fin de la partie boutons et IHM **** */

/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}
