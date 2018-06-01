#include <limits>
#include <map>
#include <sstream>
#include <vector>

#include <vtkArrowSource.h>
#include <vtkAxesActor.h>
#include <vtkCellArray.h>
#include <vtkColorTransferFunction.h>
#include <vtkGlyph3D.h>
#include <vtkImageData.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceMapper.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNamedColors.h>
#include <vtkOBJReader.h>
#include <vtkPerlinNoise.h>
#include <vtkPlaneSource.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPropPicker.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkSmartPointer.h>
#include <vtkStreamTracer.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkVersion.h>
#include <vtkXMLPolyDataWriter.h>

std::map<vtkSmartPointer<vtkActor>, std::string> mapaActores;

bool objEnEscena = true;
vtkSmartPointer<vtkActor> actorObj;
bool vectoresEnEscena = true;
std::vector<vtkSmartPointer<vtkActor>> vectores;
bool curvasEnEscena = true;
vtkSmartPointer<vtkActor> streamLineActor;

bool isInBounds(const double point[3], const double bounds[3][2]) {
  for (size_t i = 0; i < 3; i++)
    if (!(bounds[i][0] < point[i] && point[i] < bounds[i][1]))
      return false;
  return true;
}

/***********************************************
***********************************************/
class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera {
public:
  static KeyPressInteractorStyle *New();
  vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);

  virtual void OnLeftButtonDown() override {
    int *clickPos = this->GetInteractor()->GetEventPosition();
    vtkSmartPointer<vtkPropPicker> picker =
        vtkSmartPointer<vtkPropPicker>::New();
    // picker->Pick( clickPos[0], clickPos[1], 0,
    // this->GetInteractor()->GetRenderWindow()->GetRenderers()->GetFirstRenderer()
    // );
    picker->Pick(clickPos[0], clickPos[1], 0,
                 this->GetInteractor()
                     ->GetRenderWindow()
                     ->GetRenderers()
                     ->GetFirstRenderer());
    // picker->Pick( 500, 400, 0, renderer );
    // std::cout << "Picked actor static: " << picker->GetActor() << std::endl;
    vtkSmartPointer<vtkActor> seleccionado = picker->GetActor();
    if (seleccionado != 0 &&
        mapaActores.find(seleccionado) != mapaActores.end()) {
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << "====================================================="
                << std::endl;
      std::cout << "Vector: " << std::endl;
      std::cout << "-----------------------------------------------------"
                << std::endl;
      std::cout << mapaActores[picker->GetActor()];
      std::cout << "====================================================="
                << std::endl;
    }
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
  }
  virtual void OnKeyPress() override {
    // Get the keypress
    vtkRenderWindowInteractor *rwi = this->Interactor;
    std::string key = rwi->GetKeySym();

    // Handle a "normal" key
    if (key == "z") {
      if (vectoresEnEscena) {
        for (size_t i = 0; i < vectores.size(); i++) {
          vectores[i]->GetProperty()->SetOpacity(0.1);
        }
      } else {
        for (size_t i = 0; i < vectores.size(); i++) {
          vectores[i]->GetProperty()->SetOpacity(1.0);
        }
      }
      vectoresEnEscena = !vectoresEnEscena;
    } else if (key == "x") {
      if (curvasEnEscena) {
        streamLineActor->GetProperty()->SetOpacity(0.1);
      } else {
        streamLineActor->GetProperty()->SetOpacity(1.0);
      }
      curvasEnEscena = !curvasEnEscena;
    } else if (key == "c") {
      if (objEnEscena) {
        actorObj->GetProperty()->SetOpacity(0.1);
      } else {
        actorObj->GetProperty()->SetOpacity(1.0);
      }
      objEnEscena = !objEnEscena;
    }

    rwi->Render();

    // Forward events
    // vtkInteractorStyleTrackballCamera::OnKeyPress();
  }
};
vtkStandardNewMacro(KeyPressInteractorStyle);

vtkSmartPointer<vtkActor> actorVector(double x1, double y1, double z1,
                                      double x2, double y2, double z2,
                                      double minNorm, double maxNorm) {
  vtkSmartPointer<vtkArrowSource> arrowSource =
      vtkSmartPointer<vtkArrowSource>::New();
  arrowSource->SetTipResolution(20);
  arrowSource->SetShaftResolution(10);

  // Generate a random start and end point
  double startPoint[3], endPoint[3];
  vtkMath::RandomSeed(time(NULL));
  startPoint[0] = x1;
  startPoint[1] = y1;
  startPoint[2] = z1;
  endPoint[0] = x2;
  endPoint[1] = y2;
  endPoint[2] = z2;

  // Compute a basis
  double normalizedX[3];
  double normalizedY[3];
  double normalizedZ[3];

  // The X axis is a vector from start to end
  vtkMath::Subtract(endPoint, startPoint, normalizedX);
  double length = vtkMath::Norm(normalizedX);
  vtkMath::Normalize(normalizedX);

  // The Z axis is an arbitrary vector cross X
  double arbitrary[3];
  arbitrary[0] = vtkMath::Random(-10, 10);
  arbitrary[1] = vtkMath::Random(-10, 10);
  arbitrary[2] = vtkMath::Random(-10, 10);
  vtkMath::Cross(normalizedX, arbitrary, normalizedZ);
  vtkMath::Normalize(normalizedZ);

  // The Y axis is Z cross X
  vtkMath::Cross(normalizedZ, normalizedX, normalizedY);
  vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();

  // Create the direction cosine matrix
  matrix->Identity();
  for (unsigned int i = 0; i < 3; i++) {
    matrix->SetElement(i, 0, normalizedX[i]);
    matrix->SetElement(i, 1, normalizedY[i]);
    matrix->SetElement(i, 2, normalizedZ[i]);
  }

  // Apply the transforms
  vtkSmartPointer<vtkTransform> transform =
      vtkSmartPointer<vtkTransform>::New();
  transform->Translate(startPoint);
  transform->Concatenate(matrix);
  transform->Scale(length, length, length);

  // Transform the polydata
  vtkSmartPointer<vtkTransformPolyDataFilter> transformPD =
      vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformPD->SetTransform(transform);
  transformPD->SetInputConnection(arrowSource->GetOutputPort());

  // Create a mapper and actor for the arrow
  vtkSmartPointer<vtkPolyDataMapper> mapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  mapper->SetInputConnection(transformPD->GetOutputPort());
  mapper->Update();
  mapper->StaticOn();
  actor->SetMapper(mapper);
  double offset[3] = {1, 1, 1};
  double color[3] = {1, 1, 1};
  // vtkMath::Add(normalizedX, offset, color);
  vtkNew<vtkColorTransferFunction> colorTransferFunction;
  colorTransferFunction->AddRGBPoint(0, 0, 0, 0);
  colorTransferFunction->AddRGBPoint(0.1111111111111111, 0.2718098958333333, 0,
                                     0);
  colorTransferFunction->AddRGBPoint(0.2222222222222222, 0.44694010416666674, 0,
                                     0);
  colorTransferFunction->AddRGBPoint(0.3333333333333333, 0.6357421875, 0, 0);
  colorTransferFunction->AddRGBPoint(0.4444444444444444, 0.8352864583333335, 0,
                                     0);
  colorTransferFunction->AddRGBPoint(0.5555555555555556, 1, 0.19222005208333343,
                                     0);
  colorTransferFunction->AddRGBPoint(0.6666666666666666, 1, 0.49414062500000006,
                                     0);
  colorTransferFunction->AddRGBPoint(0.7777777777777778, 1, 0.6953125000000003,
                                     0);
  colorTransferFunction->AddRGBPoint(0.8888888888888888, 1, 0.8736979166666667,
                                     0);
  colorTransferFunction->AddRGBPoint(1, 1, 1, 1);
  colorTransferFunction->Build();
  double l = (length - minNorm) / (maxNorm - minNorm);
  colorTransferFunction->GetColor(l, color);
  actor->GetProperty()->SetColor(color);
  actor->GetProperty()->SetAmbientColor(color);
  actor->GetProperty()->SetAmbient(0.2);
  return actor;
}

vtkSmartPointer<vtkImageData> readImage(std::string filename,
                                        const double bounds[3][2]) {
  // Create an image
  vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();

  // Specify the size of the image data
  image->SetDimensions(10, 10, 10);

  image->AllocateScalars(VTK_FLOAT, 3);

  int *dims = image->GetDimensions();

  vtkSmartPointer<vtkPerlinNoise> genX = vtkSmartPointer<vtkPerlinNoise>::New();
  genX->SetFrequency(2.0, 2.0, 1.0);
  genX->SetAmplitude(1.0);
  vtkSmartPointer<vtkPerlinNoise> genY = vtkSmartPointer<vtkPerlinNoise>::New();
  genY->SetFrequency(1.0, 2.0, 2.0);
  genY->SetAmplitude(1.0);
  vtkSmartPointer<vtkPerlinNoise> genZ = vtkSmartPointer<vtkPerlinNoise>::New();
  genZ->SetFrequency(2.0, 1.0, 2.0);
  genZ->SetAmplitude(1.0);

  for (auto z = 0; z < dims[2]; ++z) {
    for (auto y = 0; y < dims[1]; ++y) {
      for (auto x = 0; x < dims[0]; ++x) {
        float *pixel = static_cast<float *>(image->GetScalarPointer(x, y, z));
        double point[3];
        point[0] = double(x);
        point[1] = double(y);
        point[2] = double(z);
        if (isInBounds(point, bounds)) {
          pixel[0] = 0;
          pixel[1] = 0;
          pixel[2] = 0;
        } else {
          pixel[0] = genX->EvaluateFunction(x, y, z);
          pixel[1] = genY->EvaluateFunction(x, y, z);
          pixel[2] = genZ->EvaluateFunction(x, y, z);
        }

        // for (size_t i = 0; i < 3; i++) {
        //   std::cout << "\t" << pixel[i];
        // }
        // std::cout << std::endl;
      }
    }
  }
  return image;
}

int main(int argc, char *argv[]) {
  double bounds[3][2] = {{1, 8}, {0, 4}, {1, 7}};
  std::string obj;
  if (argc == 8) {
    obj = argv[1];
    for (size_t i = 0; i < 3; i++) {
      bounds[i][0] = atof(argv[2 * i + 2]);
      bounds[i][1] = atof(argv[2 * i + 3]);
      std::cout << 2 * i + 2 << ": " << bounds[i][0] << std::endl;
      std::cout << 2 * i + 3 << ": " << bounds[i][1] << std::endl;
    }
  } else if (argc == 2) {
    obj = argv[1];
    for (size_t i = 0; i < 3; i++) {
      std::cout << 2 * i + 2 << ": " << bounds[i][0] << std::endl;
      std::cout << 2 * i + 3 << ": " << bounds[i][1] << std::endl;
    }
  } else {
    std::cout << "Error. Usage " << argv[0] << " <obj Filename>" << std::endl;
    return -1;
  }

  // ReadOBJ
  vtkSmartPointer<vtkOBJReader> readerObj =
      vtkSmartPointer<vtkOBJReader>::New();
  readerObj->SetFileName(obj.c_str());
  readerObj->Update();

  vtkSmartPointer<vtkTransform> transformObj =
      vtkSmartPointer<vtkTransform>::New();
  transformObj->Identity();
  transformObj->Translate(3, -1, 4);
  transformObj->RotateY(12);
  transformObj->RotateX(-90);
  transformObj->Scale(0.03, 0.03, 0.03);

  vtkSmartPointer<vtkTransformFilter> filterObj =
      vtkSmartPointer<vtkTransformFilter>::New();
  filterObj->SetInputConnection(readerObj->GetOutputPort());
  filterObj->SetTransform(transformObj);

  vtkSmartPointer<vtkPolyDataMapper> mapperObj =
      vtkSmartPointer<vtkPolyDataMapper>::New();
  mapperObj->SetInputConnection(filterObj->GetOutputPort());

  actorObj = vtkSmartPointer<vtkActor>::New();
  actorObj->SetMapper(mapperObj);
  // Fin ReadOBJ

  // Create an image
  vtkSmartPointer<vtkImageData> image = readImage("", bounds);
  double minNorm = std::numeric_limits<double>::infinity();
  double maxNorm = -std::numeric_limits<double>::infinity();
  int *dims = image->GetDimensions();
  for (auto z = 0; z < dims[2]; ++z) {
    for (auto y = 0; y < dims[1]; ++y) {
      for (auto x = 0; x < dims[0]; ++x) {
        float *pixel = static_cast<float *>(image->GetScalarPointer(x, y, z));
        double norm = vtkMath::Norm(pixel);
        if (norm > maxNorm) {
          maxNorm = norm;
        }
        if (norm < minNorm) {
          minNorm = norm;
        }
      }
    }
  }

  // A better way to do this is (should be tested for compatibility and
  // correctness).
  // std::cout << image->GetPointData()->GetScalars()->GetName() << std::endl;
  image->GetPointData()->SetActiveVectors("ImageScalars");

  // Setup the arrows
  vtkSmartPointer<vtkArrowSource> arrowSource =
      vtkSmartPointer<vtkArrowSource>::New();
  arrowSource->SetTipRadius(0.1);
  arrowSource->SetShaftRadius(0.025);
  arrowSource->Update();

  vtkNew<vtkMultiBlockDataSet> dataSets;
  dataSets->SetNumberOfBlocks(1);
  dataSets->SetBlock(0, image);

  vtkSmartPointer<vtkNamedColors> namedColors =
      vtkSmartPointer<vtkNamedColors>::New();

  // Source of the streamlines
  vtkSmartPointer<vtkPlaneSource> seeds =
      vtkSmartPointer<vtkPlaneSource>::New();
  seeds->SetXResolution(100);
  seeds->SetYResolution(100);
  seeds->SetOrigin(0, 0, 0);
  seeds->SetPoint1(0, 0, 9);
  seeds->SetPoint2(0, 9, 0);
  vtkSmartPointer<vtkPolyDataMapper> planeMapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
  planeMapper->SetInputConnection(seeds->GetOutputPort());
  vtkSmartPointer<vtkActor> planeActor = vtkSmartPointer<vtkActor>::New();
  planeActor->SetMapper(planeMapper);
  planeActor->GetProperty()->SetRepresentationToWireframe();
  planeActor->GetProperty()->SetColor(1, 1, 1);
  planeActor->GetProperty()->SetSpecular(0.0);
  planeActor->GetProperty()->SetDiffuse(0.0);
  planeActor->GetProperty()->SetAmbientColor(1, 1, 1);
  planeActor->GetProperty()->SetAmbient(1.0);

  // Streamline itself
  vtkSmartPointer<vtkStreamTracer> streamLine =
      vtkSmartPointer<vtkStreamTracer>::New();
  streamLine->SetInputData(dataSets);
  streamLine->SetSourceConnection(seeds->GetOutputPort());
  streamLine->SetMaximumPropagation(1000);
  streamLine->SetMaximumNumberOfSteps(1000);
  streamLine->SetInitialIntegrationStep(0.1);
  streamLine->SetIntegrationDirectionToBoth();

  // streamLine->SetStartPosition(2,-2,30);
  // as alternative to the SetSource(), which can handle multiple
  // streamlines, you can set a SINGLE streamline from
  // SetStartPosition()

  vtkSmartPointer<vtkPolyDataMapper> streamLineMapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
  streamLineMapper->SetInputConnection(streamLine->GetOutputPort());

  streamLineActor = vtkSmartPointer<vtkActor>::New();
  streamLineActor->SetMapper(streamLineMapper);
  streamLineActor->VisibilityOn();

  vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
  axes->SetTotalLength(10.0, 10.0, 10.0);

  // Setup renderer
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  for (auto z = 0; z < dims[2]; ++z) {
    for (auto y = 0; y < dims[1]; ++y) {
      for (auto x = 0; x < dims[0]; ++x) {
        float *pixel = static_cast<float *>(image->GetScalarPointer(x, y, z));
        vtkSmartPointer<vtkActor> actor = actorVector(
            double(x), double(y), double(z), double(x) + pixel[0],
            double(y) + pixel[1], double(z) + pixel[2], minNorm, maxNorm);
        vectores.push_back(actor);
        renderer->AddActor(actor);
        std::stringstream ss;
        ss << "Origen: " << x << ", " << y << ", " << z << std::endl;
        ss << "Direccion: " << pixel[0] << ", " << pixel[1] << ", " << pixel[2]
           << std::endl;
        mapaActores[actor] = ss.str();
      }
    }
  }
  renderer->AddActor(streamLineActor);
  renderer->AddActor(planeActor);
  renderer->AddActor(axes);
  renderer->AddActor(actorObj);
  renderer->ResetCamera();

  // Setup render window
  vtkSmartPointer<vtkRenderWindow> renderWindow =
      vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  // Setup render window interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
  vtkSmartPointer<KeyPressInteractorStyle> style =
      vtkSmartPointer<KeyPressInteractorStyle>::New();

  renderWindowInteractor->SetInteractorStyle(style);
  renderWindow->SetSize(1000, 1000);
  // Render and start interaction
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->Initialize();

  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
