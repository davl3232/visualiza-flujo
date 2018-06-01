#include <vtkArrowSource.h>
#include <vtkAxesActor.h>
#include <vtkCellArray.h>
#include <vtkGlyph3D.h>
#include <vtkImageData.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceMapper.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNamedColors.h>
#include <vtkPerlinNoise.h>
#include <vtkPlaneSource.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkStreamTracer.h>
#include <vtkVersion.h>
#include <vtkXMLPolyDataWriter.h>
bool isInBounds(const double point[3], const double bounds[3][2]) {
  for (size_t i = 0; i < 3; i++)
    if (!(bounds[i][0] < point[i] && point[i] < bounds[i][1]))
      return false;
  return true;
}

vtkSmartPointer<vtkImageData> readImage(std::string filename,
                                        const double bounds[3][2]) {
  // Create an image
  vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();

  // Specify the size of the image data
  image->SetDimensions(10, 10, 10);

#if VTK_MAJOR_VERSION <= 5
  image->SetNumberOfScalarComponents(3);
  image->SetScalarTypeToFloat();
  image->AllocateScalars();
#else
  image->AllocateScalars(VTK_FLOAT, 3);
#endif

  int *dims = image->GetDimensions();

  vtkSmartPointer<vtkPerlinNoise> genX = vtkSmartPointer<vtkPerlinNoise>::New();
  genX->SetFrequency(1.0, 2.0, 1.0);
  genX->SetAmplitude(1.0);
  vtkSmartPointer<vtkPerlinNoise> genY = vtkSmartPointer<vtkPerlinNoise>::New();
  genY->SetFrequency(2.0, 1.0, 1.0);
  genY->SetAmplitude(1.0);
  vtkSmartPointer<vtkPerlinNoise> genZ = vtkSmartPointer<vtkPerlinNoise>::New();
  genZ->SetFrequency(1.0, 1.0, 2.0);
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

        for (size_t i = 0; i < 3; i++) {
          std::cout << "\t" << pixel[i];
        }
        std::cout << std::endl;
      }
    }
  }
  return image;
}

int main(int argc, char *argv[]) {
  double bounds[3][2] = {{-0.5, 0.5}, {-0.5, 0.5}, {-0.5, 0.5}};
  if (argc == 7) {
    for (size_t i = 0; i < 3; i++) {
      std::cout << 2 * i + 1 << ": " << atof(argv[2 * i + 1]) << std::endl;
      std::cout << 2 * i + 2 << ": " << atof(argv[2 * i + 2]) << std::endl;
      bounds[i][0] = atof(argv[2 * i + 1]);
      bounds[i][1] = atof(argv[2 * i + 2]);
    }
  }

  // Create an image
  vtkSmartPointer<vtkImageData> image = readImage("", bounds);

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

  vtkSmartPointer<vtkGlyph3D> glyphFilter = vtkSmartPointer<vtkGlyph3D>::New();
  glyphFilter->SetSourceConnection(arrowSource->GetOutputPort());
  // glyphFilter->OrientOn();
  glyphFilter->SetScaleModeToScaleByVector();
  glyphFilter->ScalingOn();
  glyphFilter->OrientOn();
#if VTK_MAJOR_VERSION <= 5
  glyphFilter->SetInputConnection(image->GetProducerPort());
#else
  glyphFilter->SetInputData(image);
#endif
  glyphFilter->Update();

  // Create actors
  //   vtkSmartPointer<vtkImageSliceMapper> imageMapper =
  //       vtkSmartPointer<vtkImageSliceMapper>::New();
  // #if VTK_MAJOR_VERSION <= 5
  //   imageMapper->SetInputConnection(image->GetProducerPort());
  // #else
  //   imageMapper->SetInputData(image);
  // #endif

  // vtkSmartPointer<vtkImageSlice> imageSlice =
  //     vtkSmartPointer<vtkImageSlice>::New();
  // imageSlice->SetMapper(imageMapper);

  vtkSmartPointer<vtkPolyDataMapper> vectorMapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
  vectorMapper->SetInputConnection(glyphFilter->GetOutputPort());
  vtkSmartPointer<vtkActor> vectorActor = vtkSmartPointer<vtkActor>::New();
  vectorActor->SetMapper(vectorMapper);

  vtkNew<vtkMultiBlockDataSet> dataSets;
  dataSets->SetNumberOfBlocks(1);
  dataSets->SetBlock(0, image);

  vtkSmartPointer<vtkNamedColors> namedColors =
      vtkSmartPointer<vtkNamedColors>::New();

  // Source of the streamlines
  vtkSmartPointer<vtkPlaneSource> seeds =
      vtkSmartPointer<vtkPlaneSource>::New();
  seeds->SetXResolution(10);
  seeds->SetYResolution(10);
  seeds->SetOrigin(0, 0, 0);
  seeds->SetPoint1(1, 1, 9);
  seeds->SetPoint2(1, 9, 1);
  vtkSmartPointer<vtkPolyDataMapper> planeMapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
  planeMapper->SetInputConnection(seeds->GetOutputPort());
  vtkSmartPointer<vtkActor> planeActor = vtkSmartPointer<vtkActor>::New();
  planeActor->SetMapper(planeMapper);

  // Streamline itself
  vtkSmartPointer<vtkStreamTracer> streamLine =
      vtkSmartPointer<vtkStreamTracer>::New();
  streamLine->SetInputData(dataSets);
  streamLine->SetSourceConnection(seeds->GetOutputPort());
  streamLine->SetMaximumPropagation(100);
  streamLine->SetInitialIntegrationStep(0.1);
  streamLine->SetIntegrationDirectionToBoth();

  // streamLine->SetStartPosition(2,-2,30);
  // as alternative to the SetSource(), which can handle multiple
  // streamlines, you can set a SINGLE streamline from
  // SetStartPosition()

  vtkSmartPointer<vtkPolyDataMapper> streamLineMapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
  streamLineMapper->SetInputConnection(streamLine->GetOutputPort());

  vtkSmartPointer<vtkActor> streamLineActor = vtkSmartPointer<vtkActor>::New();
  streamLineActor->SetMapper(streamLineMapper);
  streamLineActor->VisibilityOn();

  vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
  axes->SetTotalLength(10.0, 10.0, 10.0);

  // Setup renderer
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  // renderer->AddViewProp(imageSlice);
  renderer->AddViewProp(vectorActor);
  renderer->AddActor(streamLineActor);
  renderer->AddActor(planeActor);
  renderer->AddActor(axes);
  renderer->ResetCamera();

  // Setup render window
  vtkSmartPointer<vtkRenderWindow> renderWindow =
      vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  // Setup render window interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
      vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle(style);
  renderWindow->SetSize(1000, 1000);
  // Render and start interaction
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->Initialize();

  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
