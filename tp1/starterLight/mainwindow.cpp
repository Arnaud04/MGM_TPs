#include "mainwindow.h"
#include "ui_mainwindow.h"

/* **** début de la partie à compléter **** */

void MainWindow::showSelections(MyMesh* _mesh)
{
    // on réinitialise les couleurs de tout le maillage
    resetAllColorsAndThickness(_mesh);

    if(vertexSelection>-1)
    {
        _mesh->set_color(_mesh->vertex_handle(vertexSelection), MyMesh::Color(0, 0, 255));
        _mesh->data(_mesh->vertex_handle(vertexSelection)).thickness = 6;
    }
    if(edgeSelection>-1)
    {
         _mesh->set_color(_mesh->edge_handle(edgeSelection), MyMesh::Color(254, 0, 0));
         _mesh->data(_mesh->edge_handle(edgeSelection)).thickness = 10;
    }
    if(faceSelection>-1)
         _mesh->set_color(_mesh->face_handle(faceSelection), MyMesh::Color(0, 255, 0));

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}


void MainWindow::showSelectionsNeighborhood(MyMesh* _mesh)
{
    // on réinitialise les couleurs de tout le maillage
    resetAllColorsAndThickness(_mesh);

    if(vertexSelection >= 0)
    {
        VertexHandle v = _mesh->vertex_handle(vertexSelection);

        _mesh->set_color(v, MyMesh::Color(255, 0, 0));

        for (MyMesh::VertexEdgeIter curentEdge = _mesh->ve_iter(v); curentEdge.is_valid(); curentEdge ++)
        {
            EdgeHandle c = *curentEdge;
            _mesh->set_color(c, MyMesh::Color(120, 0, 0));
            _mesh->data(c).thickness = 3;
        }

    }

    if(edgeSelection>=0)
    {
        EdgeHandle e = _mesh->edge_handle(edgeSelection);
        _mesh->set_color(e, MyMesh::Color(0,255,0));
        //On parcourt toutes les faces du mesh
        for(MyMesh::FaceIter f = _mesh->faces_begin(); f != _mesh->faces_end(); f++)
        {
            //On parcourt toutes les aretes sur chacune des faces
            for(MyMesh::FaceEdgeIter current_edge = _mesh->fe_iter(f); current_edge.is_valid(); current_edge ++)
            {
                EdgeHandle ce = *current_edge;
                if(e == ce) //si notre arrête selectionné égale à l'arrete courante de l'une des faces
                    _mesh->set_color(f, MyMesh::Color(0,120,0));

            }
        }
        if(faceSelection >= 0)
        {
            FaceHandle f = _mesh->face_handle(faceSelection);

            _mesh->set_color(f, MyMesh::Color(0, 0, 255));

            for (MyMesh::FaceFaceIter curentFace = _mesh->ff_iter(f); curentFace.is_valid(); curentFace ++)
            {
                FaceHandle cf = *curentFace;
                _mesh->set_color(cf, MyMesh::Color(0, 0, 120));
            }

        }
    }

    // on affiche le nouveau maillage
    displayMesh(_mesh);

}


void MainWindow::showPath(MyMesh* _mesh, int v1, int v2)
{
    // on réinitialise l'affichage
    resetAllColorsAndThickness(_mesh);

    VertexHandle start = _mesh->vertex_handle(v1);
    VertexHandle stop = _mesh->vertex_handle(v2);
    _mesh->set_color(start, MyMesh::Color(255, 195, 0));
    _mesh->set_color(stop, MyMesh::Color(255, 195, 0));
    _mesh->data(_mesh->vertex_handle(v1)).thickness = 10;
    _mesh->data(_mesh->vertex_handle(v2)).thickness = 10;

    HalfedgeHandle curr_he = mesh.halfedge_handle(start); //initialise le premier half_edge sortant

    while(mesh.from_vertex_handle(curr_he) != stop)
    {
        std::vector<int> listeSommet;
        _mesh->set_color(mesh.edge_handle(curr_he), MyMesh::Color(120, 120, 100));
        for(MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(mesh.to_vertex_handle(curr_he)); voh_it; ++voh_it) { //stocke les idx de tous les halfedge sortant du prochain sommet
            listeSommet.push_back((*voh_it).idx());
        }
        int i = rand() % listeSommet.size();
        curr_he = _mesh->halfedge_handle(listeSommet[i]);
    }

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}



void MainWindow::showBorder(MyMesh* _mesh) //<affiche bordure> met plusieurs seconde à s'executer
{
    resetAllColorsAndThickness(_mesh);

    int count = 0;
    for(MyMesh::EdgeIter currentEdge = _mesh->edges_begin(); currentEdge != _mesh->edges_end(); currentEdge ++)
    {
        EdgeHandle e = *currentEdge;
        //_mesh->set_color(e, MyMesh::Color(0,255,0));*/
        //On parcourt toutes les faces du mesh
        for(MyMesh::FaceIter f = _mesh->faces_begin(); f != _mesh->faces_end(); f++)
        {
            //On parcourt toutes les aretes sur chacune des faces
            for(MyMesh::FaceEdgeIter current_edgeFace = _mesh->fe_iter(f); current_edgeFace.is_valid(); current_edgeFace ++)
            {
                EdgeHandle ce = *current_edgeFace;
                if(e == ce) //si notre arrête selectionné égale à l'arrete courante de l'une des faces
                    count ++;
            }

        }
        if(count == 1)
            _mesh->set_color(e, MyMesh::Color(200,0,170));
            _mesh->data(edge).thickness = 6;
        count = 0;
    }
    // on affiche le nouveau maillage
    displayMesh(_mesh);
}

/* **** fin de la partie à compléter **** */


/* **** début de la partie boutons et IHM **** */

void MainWindow::on_pushButton_bordure_clicked()
{
    showBorder(&mesh);
}

void MainWindow::on_pushButton_voisinage_clicked()
{
    // changement de mode entre avec et sans voisinage
    if(modevoisinage)
    {
        ui->pushButton_voisinage->setText("Repasser en mode normal");
        modevoisinage = false;
    }
    else
    {
        ui->pushButton_voisinage->setText("Passer en mode voisinage");
        modevoisinage = true;
    }

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}


void MainWindow::on_pushButton_vertexMoins_clicked()
{
    // mise à jour de l'interface
    vertexSelection = vertexSelection - 1;
    ui->labelVertex->setText(QString::number(vertexSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_vertexPlus_clicked()
{
    // mise à jour de l'interface
    vertexSelection = vertexSelection + 1;
    ui->labelVertex->setText(QString::number(vertexSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_edgeMoins_clicked()
{
    // mise à jour de l'interface
    edgeSelection = edgeSelection - 1;
    ui->labelEdge->setText(QString::number(edgeSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_edgePlus_clicked()
{
    // mise à jour de l'interface
    edgeSelection = edgeSelection + 1;
    ui->labelEdge->setText(QString::number(edgeSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_faceMoins_clicked()
{
    // mise à jour de l'interface
    faceSelection = faceSelection - 1;
    ui->labelFace->setText(QString::number(faceSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_facePlus_clicked()
{
    // mise à jour de l'interface
    faceSelection = faceSelection + 1;
    ui->labelFace->setText(QString::number(faceSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_afficherChemin_clicked()
{
    // on récupère les sommets de départ et d'arrivée
    int indexV1 = ui->spinBox_v1_chemin->value();
    int indexV2 = ui->spinBox_v2_chemin->value();

    showPath(&mesh, indexV1, indexV2);
}


void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

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
void MainWindow::displayMesh(MyMesh* _mesh)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
    MyMesh::ConstFaceVertexIter fvIt;
    int i = 0;
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

