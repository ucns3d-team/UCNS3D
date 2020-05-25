// Adaptor for getting Fortran simulation code into ParaView CoProcessor.

// CoProcessor specific headers
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkDoubleArray.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkCellData.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"

// Fortran specific header
#include "vtkCPPythonAdaptorAPI.h"

// These will be called from the Fortran "glue" code"
// Completely dependent on data layout, structured vs. unstructured, etc.
// since VTK/ParaView uses different internal layouts for each.

// Creates the data container for the CoProcessor.
extern "C" void createcpimagedata_(int* nxstart, int* nxend, int* nx, int* ny, int* nz)
{
  if (!vtkCPPythonAdaptorAPI::GetCoProcessorData())
  {
    vtkGenericWarningMacro("Unable to access CoProcessorData.");
    return;
  }

  // The simulation grid is a 3-dimensional topologically and geometrically
  // regular grid. In VTK/ParaView, this is considered an image data set.
  vtkSmartPointer<vtkImageData> grid = vtkSmartPointer<vtkImageData>::New();

  grid->SetExtent(*nxstart - 1, *nxend - 1, 0, *ny - 1, 0, *nz - 1);

  // Name should be consistent between here, Fortran and Python client script.
  vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input")->SetGrid(grid);
  vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input")->SetWholeExtent(
    0, *nx - 1, 0, *ny - 1, 0, *nz - 1);
}

extern "C" void testfunction_(double pointSet[][3], int* pointSetSize, int vtkCellId[][8], int* vtkCellIdSetSize, int* rank, int* numprocs)
{
  // if (*rank == 0){
  if (!vtkCPPythonAdaptorAPI::GetCoProcessorData())
  {
    vtkGenericWarningMacro("Unable to access CoProcessorData.");
    return;
  }

  vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
  cellIdList->SetNumberOfIds(8);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(*pointSetSize);

  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();

  static vtkIdType ids[8];
  double point[3];

  std::cout << "Process " << *rank << " with " <<  *numprocs << " processes involved."<< "\n";
  std::cout << "Checking size of pointSet: " << *pointSetSize << "\n";
  std::cout << "Checking size of vtkCellId: " << *vtkCellIdSetSize << "\n\n";

  // Testing pointSet output
  // if (*rank == 0){
  //   for(int i = 0; i < *pointSetSize; i++){
  //     std::cout << i+1 << ' ' << pointSet[i][0] << std::setprecision(10) << ' '
  //               << pointSet[i][1] << std::setprecision(10) << ' '
  //               << pointSet[i][2] << std::setprecision(10) << '\n';
  //   }
  // }
  // Testing vtkCellId output
  // if(*rank == 0){
    // for(int j = 0; j < *vtkCellIdSetSize; j++){
    //   for(int k = 0; k < 8; k++){
    //     connectivity_array[j*8 + k] = vtkCellId[j][k];
    //     if(k == 0){
    //       offsets_array[j] = j*8;
    //     }
    //   }
    // }

  // Serial process of registering point set
  for(int i=0; i < *pointSetSize; i++){
    points->SetPoint(i, pointSet[i][0], pointSet[i][1], pointSet[i][2]);
  }

  unstructuredGrid->Allocate(*vtkCellIdSetSize);
  for(int k = 0; k < *vtkCellIdSetSize; k++){
    for(int i=0; i < 8; i++){
      ids[i] = vtkCellId[k][i];
    }
    unstructuredGrid->InsertNextCell(12, 8, ids);
  }

  // if (*rank == 0){
    // for(int i = 0; i < *pointSetSize; i++){
      std::cout << "Total number of cells : " << unstructuredGrid->GetNumberOfCells() << "\n";
    // }
  // }

  unstructuredGrid->SetPoints(points);
  //
  vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input")->SetGrid(unstructuredGrid);
}

// Add field(s) to the data container.
// Separate from above because this will be dynamic, grid is static.
// By hand name mangling for fortran.
extern "C" void addfield_(double* scalars, char* name)
{
  vtkCPInputDataDescription* idd =
    vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input");

  // vtkImageData* image = vtkImageData::SafeDownCast(idd->GetGrid());

  vtkUnstructuredGrid* unstructuredGrid = vtkUnstructuredGrid::SafeDownCast(idd->GetGrid());

  if (!unstructuredGrid)
  {
    vtkGenericWarningMacro("No adaptor grid to attach field data to.");
    return;
  }

  // field name must match that in the fortran code.
  if (idd->IsFieldNeeded(name, vtkDataObject::CELL))
  {
    vtkSmartPointer<vtkDoubleArray> field = vtkSmartPointer<vtkDoubleArray>::New();
    field->SetNumberOfComponents(1);
    field->SetName(name);
    field->SetArray(scalars, unstructuredGrid->GetNumberOfCells(), 1);
    unstructuredGrid->GetCellData()->AddArray(field);
  }
}
