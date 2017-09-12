#define vtkRenderingCore_AUTOINIT 3(vtkInteractionStyle, vtkRenderingFreeType, vtkRenderingOpenGL2)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL2)

#include <vtkSmartPointer.h>
#include <vtkImageViewer2.h>
#include <vtkImageActor.h>
#include <vtkPNGReader.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkImageData.h>
#include <iostream>
#include "Vector3.h"


/*
	Calculates a new point using bilinear interpolation
	@param P1
	@param P2
	@param P3
	@param P4
	@param u
	@param v
	@return
*/
Vector3* bilinearInterpolation(Vector3* P1, Vector3* P2, Vector3* P3, Vector3* P4, double u, double v) 
{
	// Calculate the middle vectors
	Vector3* t1 = P2->sub(P1);
	Vector3* t2 = P4->sub(P1);
	Vector3* t3 = P1->sub(P2)->add(P3)->sub(P4);

	// Calculate the location of the point
	return P1->add(t1->multiply(u))->add(t2->multiply(v))->add(t3->multiply(u * v));
}


/*
	Calculates a new point using barycentric coordinates
	@param P1
	@param P2
	@param alpha
	@param beta
	@param gamma
	@return
*/
Vector3* barycentricCoordinates(Vector3* P1, Vector3* P2, Vector3* P3, double alpha, double beta, double gamma)
{
	// Calculate the location of the point
	return P1->multiply(alpha)->add(P2->multiply(beta))->add(P3->multiply(gamma));
}


/*
	Returns an interpolated version of the image by using bilinear interpolation in nxn regions
	@param imageData
	@param n
	@return
*/
vtkSmartPointer<vtkImageData> useBilinearInterpolation(vtkSmartPointer<vtkImageData> imageData, int n)
{
	// Get the dimensions of the image
	int* dims;
	dims = imageData->GetDimensions();

	// Pointer to the new image where interpolated values will be drawn
	// NOTE: The scalar type of the original image is unsigned char with 3 components
	vtkSmartPointer<vtkImageData> newImage = vtkSmartPointer<vtkImageData>::New();
	newImage->SetDimensions(dims[0], dims[1], dims[2]);
	newImage->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

	// Get the size of each cell length
	int size = dims[0] / n;

	// Initialize the sum of squared difference
	double difference = 0.0;

	// Traverse the image by squared regions
	for (int y = 0; y < dims[1]; y += size)
	{
		for (int x = 0; x < dims[0]; x += size)
		{
			// Get the (x,y) location for the next cells
			int next_cell_x = (x + size - 1 >= dims[0]) ? dims[0] - 1 : x + size - 1;
			int next_cell_y = (y + size - 1 >= dims[1]) ? dims[1] - 1 : y + size - 1;

			// Define the color for each square vertex
			unsigned char* pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(x, y, 0));
			Vector3* C1 = new Vector3((int)pixel[0], (int)pixel[1], (int)pixel[2]);

			pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(next_cell_x, y, 0));
			Vector3* C2 = new Vector3((int)pixel[0], (int)pixel[1], (int)pixel[2]);

			pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(next_cell_x, next_cell_y, 0));
			Vector3* C3 = new Vector3((int)pixel[0], (int)pixel[1], (int)pixel[2]);

			pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(x, next_cell_y, 0));
			Vector3* C4 = new Vector3((int)pixel[0], (int)pixel[1], (int)pixel[2]);

			// Traverse each squared regions and calculate the middle points color
			for (int py = 0; py < size; py += 1)
			{
				for (int px = 0; px < size; px += 1)
				{
					// Get the location of the current pixel
					int pixel_x = (x + px >= dims[0]) ? dims[0] - 1 : x + px;
					int pixel_y = (y + py >= dims[1]) ? dims[1] - 1 : y + py;

					// Get the original color of the current pixel
					pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(pixel_x, pixel_y, 0));
					Vector3* originalColor = new Vector3((double)pixel[0], (double)pixel[1], (double)pixel[2]);

					// Calculate the (u, v) values for the current pixel
					double u = double(px) / double(size);
					double v = double(py) / double(size);

					// Calculate the color of the point in the squared region using bilinear interpolation
					Vector3* newColor = bilinearInterpolation(C1, C2, C3, C4, u, v);

					// Set the new RGB values
					pixel = static_cast<unsigned char*>(newImage->GetScalarPointer(pixel_x, pixel_y, 0));
					pixel[0] = (unsigned char)((int)newColor->x);
					pixel[1] = (unsigned char)((int)newColor->y);
					pixel[2] = (unsigned char)((int)newColor->z);
					Vector3* C = new Vector3((double)pixel[0], (double)pixel[1], (double)pixel[2]);

					// Accumulate the squared distance
					difference += originalColor->squareDistance(C);
				}
			}
		}
	}

	// Write the error
	std::cout << difference << std::endl;

	// Return the pointer to the new image
	return newImage;
}


/*
	Returns an interpolated version of the image by using barycentric coordinates in nxn regions
	@param imageData
	@param n
	@return
*/
vtkSmartPointer<vtkImageData> useBarycentricCoordinates(vtkSmartPointer<vtkImageData> imageData, int n)
{
	// Get the dimensions of the image
	int* dims;
	dims = imageData->GetDimensions();

	// Pointer to the new image where interpolated values will be drawn
	// NOTE: The scalar type of the original image is unsigned char with 3 components
	vtkSmartPointer<vtkImageData> newImage = vtkSmartPointer<vtkImageData>::New();
	newImage->SetDimensions(dims[0], dims[1], dims[2]);
	newImage->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

	// Get the size of each cell length
	int size = dims[0] / n;

	// Initialize the sum of squared difference
	double difference = 0.0;

	// Traverse the image by squared regions
	for (int y = 0; y < dims[1]; y += size)
	{
		for (int x = 0; x < dims[0]; x += size)
		{
			// Get the (x,y) location for the next cells
			int next_cell_x = (x + size - 1 >= dims[0]) ? dims[0] - 1 : x + size - 1;
			int next_cell_y = (y + size - 1 >= dims[1]) ? dims[1] - 1 : y + size - 1;

			// Define the color for each square vertex
			Vector3* P1 = new Vector3(x, y, 0);
			unsigned char* pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(x, y, 0));
			Vector3* C1 = new Vector3((int)pixel[0], (int)pixel[1], (int)pixel[2]);

			Vector3* P2 = new Vector3(x + size - 1, y, 0);
			pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(next_cell_x, y, 0));
			Vector3* C2 = new Vector3((int)pixel[0], (int)pixel[1], (int)pixel[2]);

			Vector3* P3 = new Vector3(x + size - 1, y + size - 1, 0);
			pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(next_cell_x, next_cell_y, 0));
			Vector3* C3 = new Vector3((int)pixel[0], (int)pixel[1], (int)pixel[2]);

			Vector3* P4 = new Vector3(x, y + size - 1, 0);
			pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(x, next_cell_y, 0));
			Vector3* C4 = new Vector3((int)pixel[0], (int)pixel[1], (int)pixel[2]);

			// Traverse each squared regions and calculate the middle points color
			for (int py = 0; py < size; py += 1)
			{
				for (int px = 0; px < size; px += 1)
				{
					// Get the location of the current pixel
					int pixel_x = (x + px >= dims[0]) ? dims[0] - 1 : x + px;
					int pixel_y = (y + py >= dims[1]) ? dims[1] - 1 : y + py;

					// Get the original color of the current pixel
					pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(pixel_x, pixel_y, 0));
					Vector3* originalColor = new Vector3((double)pixel[0], (double)pixel[1], (double)pixel[2]);

					// Define the location point for the current pixel
					Vector3* P = new Vector3(pixel_x, pixel_y, 0);

					// The pointer to the vector with the color for the current pixel
					Vector3* newColor = NULL;

					// if px > py then pixel is in the lower triangle of the square; otherwise it is in the 
					// upper triangle of the square
					if (px >= py) 
					{
						// Find the barycentric coordinates for the current pixel
						Vector3* v0 = P2->sub(P1);
						Vector3* v1 = P3->sub(P1);
						Vector3* v2 = P->sub(P1);
						double d00 = v0->dot(v0);
						double d01 = v0->dot(v1);
						double d11 = v1->dot(v1);
						double d20 = v2->dot(v0);
						double d21 = v2->dot(v1);
						double denom = d00 * d11 - d01 * d01;
						double b = (d11 * d20 - d01 * d21) / denom;
						double g = (d00 * d21 - d01 * d20) / denom;
						double a = 1.0 - b - g;

						// Calculate the color of the point in the squared region using bilinear interpolation
						newColor = barycentricCoordinates(C1, C2, C3, a, b, g);
					}
					else 
					{
						// Find the barycentric coordinates for the current pixel
						Vector3* v0 = P4->sub(P1);
						Vector3* v1 = P3->sub(P1);
						Vector3* v2 = P->sub(P1);
						double d00 = v0->dot(v0);
						double d01 = v0->dot(v1);
						double d11 = v1->dot(v1);
						double d20 = v2->dot(v0);
						double d21 = v2->dot(v1);
						double denom = d00 * d11 - d01 * d01;
						double b = (d11 * d20 - d01 * d21) / denom;
						double g = (d00 * d21 - d01 * d20) / denom;
						double a = 1.0 - b - g;

						// Calculate the color of the point in the squared region using bilinear interpolation
						newColor = barycentricCoordinates(C1, C4, C3, a, g, b);
					}

					// Set the new RGB values
					pixel = static_cast<unsigned char*>(newImage->GetScalarPointer(pixel_x, pixel_y, 0));
					pixel[0] = (unsigned char)((int)newColor->x);
					pixel[1] = (unsigned char)((int)newColor->y);
					pixel[2] = (unsigned char)((int)newColor->z);
					Vector3* C = new Vector3((double)pixel[0], (double)pixel[1], (double)pixel[2]);

					// Accumulate the squared distance
					difference += originalColor->squareDistance(C);
				}
			}
		}
	}

	// Write the error
	std::cout << difference << std::endl;

	// Return the pointer to the new image
	return newImage;
}


/*
	Returns an interpolated version of the image by using Shepard in nxn regions
	@param imageData
	@param n
	@return
*/
vtkSmartPointer<vtkImageData> useShepard(vtkSmartPointer<vtkImageData> imageData, int n)
{
	// Get the dimensions of the image
	int* dims;
	dims = imageData->GetDimensions();

	// Pointer to the new image where interpolated values will be drawn
	// NOTE: The scalar type of the original image is unsigned char with 3 components
	vtkSmartPointer<vtkImageData> newImage = vtkSmartPointer<vtkImageData>::New();
	newImage->SetDimensions(dims[0], dims[1], dims[2]);
	newImage->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

	// Get the size of each cell length
	int size = dims[0] / n;

	// Initialize the sum of squared difference
	double difference = 0.0;

	// Traverse the image by squared regions
	for (int y = 0; y < dims[1]; y += size)
	{
		for (int x = 0; x < dims[0]; x += size)
		{
			// Get the (x,y) location for the next cells
			int next_cell_x = (x + size - 1 >= dims[0]) ? dims[0] - 1 : x + size - 1;
			int next_cell_y = (y + size - 1 >= dims[1]) ? dims[1] - 1 : y + size - 1;

			// Define the location and color for each square vertex
			Vector3* P1 = new Vector3(x, y, 0);
			unsigned char* pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(x, y, 0));
			Vector3* C1 = new Vector3((int)pixel[0], (int)pixel[1], (int)pixel[2]);

			Vector3* P2 = new Vector3(x + size - 1, y, 0);
			pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(next_cell_x, y, 0));
			Vector3* C2 = new Vector3((int)pixel[0], (int)pixel[1], (int)pixel[2]);

			Vector3* P3 = new Vector3(x + size - 1, y + size - 1, 0);
			pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(next_cell_x, next_cell_y, 0));
			Vector3* C3 = new Vector3((int)pixel[0], (int)pixel[1], (int)pixel[2]);

			Vector3* P4 = new Vector3(x, y + size - 1, 0);
			pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(x, next_cell_y, 0));
			Vector3* C4 = new Vector3((int)pixel[0], (int)pixel[1], (int)pixel[2]);

			// Traverse each squared regions and calculate the middle points color
			for (int py = 0; py < size; py += 1)
			{
				for (int px = 0; px < size; px += 1)
				{
					// Get the location of the current pixel
					int pixel_x = (x + px >= dims[0]) ? dims[0] - 1 : x + px;
					int pixel_y = (y + py >= dims[1]) ? dims[1] - 1 : y + py;

					// Get the original color of the current pixel
					pixel = static_cast<unsigned char*>(imageData->GetScalarPointer(pixel_x, pixel_y, 0));
					Vector3* originalColor = new Vector3((double)pixel[0], (double)pixel[1], (double)pixel[2]);

					// Define the location of the current pixel
					Vector3* P = new Vector3(pixel_x, pixel_y, 0);

					// Calculate the four weights
					double w1 = exp(-pow(P->sub(P1)->magnitude(), 2));
					double w2 = exp(-pow(P->sub(P2)->magnitude(), 2));
					double w3 = exp(-pow(P->sub(P3)->magnitude(), 2));
					double w4 = exp(-pow(P->sub(P4)->magnitude(), 2));

					// Calculate Z
					double Z = w1 + w2 + w3 + w4;

					// Calculate the color of the current pixel
					Vector3* t1 = C1->multiply(w1);
					Vector3* t2 = C2->multiply(w2);
					Vector3* t3 = C3->multiply(w3);
					Vector3* t4 = C4->multiply(w4);
					Vector3* newColor = t1->add(t2)->add(t3)->add(t4)->multiply(1.0 / Z);

					// Set the new RGB values
					pixel = static_cast<unsigned char*>(newImage->GetScalarPointer(pixel_x, pixel_y, 0));
					pixel[0] = (unsigned char)((int)newColor->x);
					pixel[1] = (unsigned char)((int)newColor->y);
					pixel[2] = (unsigned char)((int)newColor->z);
					Vector3* C = new Vector3((double)pixel[0], (double)pixel[1], (double)pixel[2]);

					// Accumulate the squared distance
					difference += originalColor->squareDistance(C);
				}
			}
		}
	}

	// Write the error
	std::cout << difference << std::endl;

	// Return the pointer to the new image
	return newImage;
}


/*
	The main function
*/
int main(int argc, char** argv)
{
	// Read the size of the grid
	// NOTE: The logic is based on cells instead of point grid. In order to make the 
	// logic point based then subtract 1 to the given n value
	// e.g., a 3x3 cell grid is a 4x4 point grid
	int n = 2;
	std::cout << "n: ";
	std::cin >> n;
	n -= 1;

	// Initialize the image reader and read the image
	vtkSmartPointer<vtkPNGReader> imageReader = vtkSmartPointer<vtkPNGReader>::New();
	imageReader->SetFileName("xid-35316670_1.png");
	imageReader->Update();

	// Get the image data from the read image. This is the pointer to the original image
	vtkSmartPointer<vtkImageData> imageData = imageReader->GetOutput();
	
	// Initialize the renderer window
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();

	// Initialize the window interactor
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	
	// Define the margins of each viewport
	double xmins[4] = { 0.0, 0.5, 0.0, 0.5 };
	double ymins[4] = { 0.5, 0.5, 0.0, 0.0 };
	double xmaxs[4] = { 0.5, 1.0, 0.5, 1.0 };
	double ymaxs[4] = { 1.0, 1.0, 0.5, 0.5 };

	// Define the 4 images: original, bilinear, barycentric and Shepard
	// NOTE: Here are obtained the interpolated images
	vtkSmartPointer<vtkImageData> images[4];
	images[0] = imageReader->GetOutput();
	images[1] = useBilinearInterpolation(imageData, n);
	images[2] = useBarycentricCoordinates(imageData, n);
	images[3] = useShepard(imageData, n);

	// Define each one of the viewports
	for (int i = 0; i < 4; i += 1) 
	{
		// Initialize the renderer
		vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

		// Add the renderer to the window
		renderWindow->AddRenderer(renderer);

		// Define the renderer viewport
		renderer->SetViewport(xmins[i], ymins[i], xmaxs[i], ymaxs[i]);

		// Initialize the input actor
		vtkSmartPointer<vtkImageActor> imageActor = vtkSmartPointer<vtkImageActor>::New();
		imageActor->SetInputData(images[i]);

		// Set the image actor to the renderer
		renderer->AddActor(imageActor);
	}

	// Set up the renderer window
	renderWindow->SetSize(600, 600);
	renderWindow->Render();

	// Start the program
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}
