#include <vtkArrowSource.h>
#include <vtkCellArray.h>
#include <vtkGlyph3D.h>
#include <vtkImageData.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceMapper.h>
#include <vtkInteractorStyleTrackballCamera.h>
// #include <vtkInteractorStyleImage.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>
#include <vtkXMLPolyDataWriter.h>

int main(int, char *[]) {
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

  // Zero the image
  for (auto z = 0; z < dims[2]; ++z) {
    for (auto y = 0; y < dims[1]; ++y) {
      for (auto x = 0; x < dims[0]; ++x) {
        float *pixel = static_cast<float *>(image->GetScalarPointer(x, y, z));
        pixel[0] = x / float(dims[0]);
        pixel[1] = y / float(dims[1]);
        pixel[2] = z / float(dims[2]);
      }
    }
  }

  // A better way to do this is (should be tested for compatibility and
  // correctness).
  // std::cout << image->GetPointData()->GetScalars()->GetName() << std::endl;
  image->GetPointData()->SetActiveVectors(
      image->GetPointData()->GetScalars()->GetName());
  // image->GetPointData()->SetActiveVectors("ImageScalars");

  // Setup the arrows
  vtkSmartPointer<vtkArrowSource> arrowSource =
      vtkSmartPointer<vtkArrowSource>::New();
  arrowSource->SetTipRadius(0.1);
  arrowSource->SetShaftRadius(0.025);
  arrowSource->Update();

  vtkSmartPointer<vtkGlyph3D> glyphFilter = vtkSmartPointer<vtkGlyph3D>::New();
  glyphFilter->SetSourceConnection(arrowSource->GetOutputPort());
  glyphFilter->OrientOn();
  glyphFilter->SetVectorModeToUseVector();
#if VTK_MAJOR_VERSION <= 5
  glyphFilter->SetInputConnection(image->GetProducerPort());
#else
  glyphFilter->SetInputData(image);
#endif
  glyphFilter->Update();

  // Create actors
  vtkSmartPointer<vtkImageSliceMapper> imageMapper =
      vtkSmartPointer<vtkImageSliceMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  imageMapper->SetInputConnection(image->GetProducerPort());
#else
  imageMapper->SetInputData(image);
#endif

  vtkSmartPointer<vtkImageSlice> imageSlice =
      vtkSmartPointer<vtkImageSlice>::New();
  imageSlice->SetMapper(imageMapper);

  vtkSmartPointer<vtkPolyDataMapper> vectorMapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
  vectorMapper->SetInputConnection(glyphFilter->GetOutputPort());
  vtkSmartPointer<vtkActor> vectorActor = vtkSmartPointer<vtkActor>::New();
  vectorActor->SetMapper(vectorMapper);

  // Setup renderer
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->AddViewProp(imageSlice);
  renderer->AddViewProp(vectorActor);
  renderer->ResetCamera();

  // Setup render window
  vtkSmartPointer<vtkRenderWindow> renderWindow =
      vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  // Setup render window interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
  // vtkSmartPointer<vtkInteractorStyleImage> style =
  //     vtkSmartPointer<vtkInteractorStyleImage>::New();
  // renderWindowInteractor->SetInteractorStyle(style);
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
      vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle(style);
  // Render and start interaction
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->Initialize();

  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
