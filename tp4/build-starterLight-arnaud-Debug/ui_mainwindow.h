/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.12.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QFrame>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPlainTextEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include <meshviewerwidget.h>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QHBoxLayout *horizontalLayout;
    QWidget *widget_2;
    QVBoxLayout *verticalLayout;
    QPushButton *pushButton_chargement;
    QFrame *frame;
    QVBoxLayout *verticalLayout_3;
    QWidget *widget_3;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_5;
    QPushButton *pushButton_edgeMoins;
    QLabel *labelEdge;
    QPushButton *pushButton_edgePlus;
    QPlainTextEdit *infoEdgeSelection;
    QPushButton *pushButton_delSelEdge;
    QFrame *frame_2;
    QVBoxLayout *verticalLayout_2;
    QWidget *widget;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label;
    QComboBox *comboBox;
    QWidget *widget_4;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label_2;
    QSlider *horizontalSlider;
    QPushButton *pushButton_decimate;
    QSpacerItem *verticalSpacer;
    MeshViewerWidget *displayWidget;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(632, 408);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        horizontalLayout = new QHBoxLayout(centralWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        widget_2 = new QWidget(centralWidget);
        widget_2->setObjectName(QString::fromUtf8("widget_2"));
        QSizePolicy sizePolicy(QSizePolicy::Maximum, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(150);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(widget_2->sizePolicy().hasHeightForWidth());
        widget_2->setSizePolicy(sizePolicy);
        widget_2->setMinimumSize(QSize(150, 0));
        verticalLayout = new QVBoxLayout(widget_2);
        verticalLayout->setSpacing(4);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(3, 3, 3, 3);
        pushButton_chargement = new QPushButton(widget_2);
        pushButton_chargement->setObjectName(QString::fromUtf8("pushButton_chargement"));
        pushButton_chargement->setMinimumSize(QSize(200, 0));

        verticalLayout->addWidget(pushButton_chargement);

        frame = new QFrame(widget_2);
        frame->setObjectName(QString::fromUtf8("frame"));
        frame->setFrameShape(QFrame::StyledPanel);
        frame->setFrameShadow(QFrame::Raised);
        verticalLayout_3 = new QVBoxLayout(frame);
        verticalLayout_3->setSpacing(1);
        verticalLayout_3->setContentsMargins(11, 11, 11, 11);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        verticalLayout_3->setContentsMargins(1, 1, 1, 1);
        widget_3 = new QWidget(frame);
        widget_3->setObjectName(QString::fromUtf8("widget_3"));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(widget_3->sizePolicy().hasHeightForWidth());
        widget_3->setSizePolicy(sizePolicy1);
        horizontalLayout_3 = new QHBoxLayout(widget_3);
        horizontalLayout_3->setSpacing(1);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(1, 1, 1, 1);
        label_5 = new QLabel(widget_3);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        QSizePolicy sizePolicy2(QSizePolicy::Maximum, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(label_5->sizePolicy().hasHeightForWidth());
        label_5->setSizePolicy(sizePolicy2);

        horizontalLayout_3->addWidget(label_5);

        pushButton_edgeMoins = new QPushButton(widget_3);
        pushButton_edgeMoins->setObjectName(QString::fromUtf8("pushButton_edgeMoins"));
        QSizePolicy sizePolicy3(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(pushButton_edgeMoins->sizePolicy().hasHeightForWidth());
        pushButton_edgeMoins->setSizePolicy(sizePolicy3);

        horizontalLayout_3->addWidget(pushButton_edgeMoins);

        labelEdge = new QLabel(widget_3);
        labelEdge->setObjectName(QString::fromUtf8("labelEdge"));
        sizePolicy1.setHeightForWidth(labelEdge->sizePolicy().hasHeightForWidth());
        labelEdge->setSizePolicy(sizePolicy1);
        labelEdge->setAlignment(Qt::AlignCenter);

        horizontalLayout_3->addWidget(labelEdge);

        pushButton_edgePlus = new QPushButton(widget_3);
        pushButton_edgePlus->setObjectName(QString::fromUtf8("pushButton_edgePlus"));
        sizePolicy3.setHeightForWidth(pushButton_edgePlus->sizePolicy().hasHeightForWidth());
        pushButton_edgePlus->setSizePolicy(sizePolicy3);
        pushButton_edgePlus->setMinimumSize(QSize(10, 0));

        horizontalLayout_3->addWidget(pushButton_edgePlus);


        verticalLayout_3->addWidget(widget_3);

        infoEdgeSelection = new QPlainTextEdit(frame);
        infoEdgeSelection->setObjectName(QString::fromUtf8("infoEdgeSelection"));
        QSizePolicy sizePolicy4(QSizePolicy::Expanding, QSizePolicy::Maximum);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(infoEdgeSelection->sizePolicy().hasHeightForWidth());
        infoEdgeSelection->setSizePolicy(sizePolicy4);
        infoEdgeSelection->setMaximumSize(QSize(200, 70));
        infoEdgeSelection->setStyleSheet(QString::fromUtf8("font: 11pt \"Lucida Console\";"));

        verticalLayout_3->addWidget(infoEdgeSelection);

        pushButton_delSelEdge = new QPushButton(frame);
        pushButton_delSelEdge->setObjectName(QString::fromUtf8("pushButton_delSelEdge"));

        verticalLayout_3->addWidget(pushButton_delSelEdge);


        verticalLayout->addWidget(frame);

        frame_2 = new QFrame(widget_2);
        frame_2->setObjectName(QString::fromUtf8("frame_2"));
        frame_2->setFrameShape(QFrame::StyledPanel);
        frame_2->setFrameShadow(QFrame::Raised);
        verticalLayout_2 = new QVBoxLayout(frame_2);
        verticalLayout_2->setSpacing(1);
        verticalLayout_2->setContentsMargins(11, 11, 11, 11);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        verticalLayout_2->setContentsMargins(1, 3, 1, 1);
        widget = new QWidget(frame_2);
        widget->setObjectName(QString::fromUtf8("widget"));
        horizontalLayout_2 = new QHBoxLayout(widget);
        horizontalLayout_2->setSpacing(1);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(1, 1, 1, 1);
        label = new QLabel(widget);
        label->setObjectName(QString::fromUtf8("label"));
        QSizePolicy sizePolicy5(QSizePolicy::Maximum, QSizePolicy::Maximum);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy5);

        horizontalLayout_2->addWidget(label);

        comboBox = new QComboBox(widget);
        comboBox->addItem(QString());
        comboBox->addItem(QString());
        comboBox->addItem(QString());
        comboBox->addItem(QString());
        comboBox->setObjectName(QString::fromUtf8("comboBox"));

        horizontalLayout_2->addWidget(comboBox);


        verticalLayout_2->addWidget(widget);

        widget_4 = new QWidget(frame_2);
        widget_4->setObjectName(QString::fromUtf8("widget_4"));
        horizontalLayout_4 = new QHBoxLayout(widget_4);
        horizontalLayout_4->setSpacing(1);
        horizontalLayout_4->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        horizontalLayout_4->setContentsMargins(1, 1, 1, 1);
        label_2 = new QLabel(widget_4);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        horizontalLayout_4->addWidget(label_2);

        horizontalSlider = new QSlider(widget_4);
        horizontalSlider->setObjectName(QString::fromUtf8("horizontalSlider"));
        horizontalSlider->setMinimum(1);
        horizontalSlider->setValue(50);
        horizontalSlider->setOrientation(Qt::Horizontal);

        horizontalLayout_4->addWidget(horizontalSlider);


        verticalLayout_2->addWidget(widget_4);

        pushButton_decimate = new QPushButton(frame_2);
        pushButton_decimate->setObjectName(QString::fromUtf8("pushButton_decimate"));

        verticalLayout_2->addWidget(pushButton_decimate);


        verticalLayout->addWidget(frame_2);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);


        horizontalLayout->addWidget(widget_2);

        displayWidget = new MeshViewerWidget(centralWidget);
        displayWidget->setObjectName(QString::fromUtf8("displayWidget"));

        horizontalLayout->addWidget(displayWidget);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 632, 22));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", nullptr));
        pushButton_chargement->setText(QApplication::translate("MainWindow", "Charger OBJ", nullptr));
        label_5->setText(QApplication::translate("MainWindow", "e : ", nullptr));
        pushButton_edgeMoins->setText(QApplication::translate("MainWindow", "-", nullptr));
        labelEdge->setText(QApplication::translate("MainWindow", "-1", nullptr));
        pushButton_edgePlus->setText(QApplication::translate("MainWindow", "+", nullptr));
        infoEdgeSelection->setPlainText(QApplication::translate("MainWindow", "Surface :\n"
"C1 : \n"
"C2 : \n"
"C3 :", nullptr));
        pushButton_delSelEdge->setText(QApplication::translate("MainWindow", "Supprimer l'ar\303\252te", nullptr));
        label->setText(QApplication::translate("MainWindow", "Crit\303\250re : ", nullptr));
        comboBox->setItemText(0, QApplication::translate("MainWindow", "Al\303\251atoire", nullptr));
        comboBox->setItemText(1, QApplication::translate("MainWindow", "Par taille", nullptr));
        comboBox->setItemText(2, QApplication::translate("MainWindow", "Par angle", nullptr));
        comboBox->setItemText(3, QApplication::translate("MainWindow", "Par plan\303\251it\303\251", nullptr));

        label_2->setText(QApplication::translate("MainWindow", "% : ", nullptr));
        pushButton_decimate->setText(QApplication::translate("MainWindow", "Simplifier avec le crit\303\250re", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
