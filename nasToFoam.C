/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    nasToFoam

Group
    grpMeshConversionUtilities

Description
    nastran dat format mesh conversion.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "IFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Dat file format.
enum struct FORMAT : label
{
    FREE = 0,
    SMALL = 8,
    LARGE = 16
};
// File format, small by default.
FORMAT format = FORMAT::SMALL;

// Last extracted entry kw
string entryBuff;
// fwd get entry
string& getEntry(IFstream& is);


bool findBulk(IFstream& is)
{
    string line;
    while (is.good())
    {
        // Skip commented lines
        if (is.peek() == '$') 
        {
            is.getLine(nullptr);
            continue;
        };

        is.getLine(line);
        if (line.starts_with("BEGIN BULK")) return true;
    }
    return false;
}

// Skip everything on this line, and next lines if we have a multiline.
void skipLine(IFstream& is)
{
    string line;
    while (true)
    {
        is.getLine(line);

        // Remove stupid carriage return char...
        line.erase(std::remove(line.begin(), line.end(), 13),
                 line.end());

        // It could be "if (!line.endsWith("+\n")) break;" but it is a little
        // confusing...
        if (line.ends_with('+')) continue;
        else break;
    }

}

string getColumn(IFstream& is, label size = label(format))
{
    label i = 0;
    char buff[63];  // More than enough, get size no chars

    if (size == 0) // Increase for free format to avoid overflow
    {
        size = 63;
    }

    while (i < size)
    {
        is.get(buff[i]);
        // Stop at new line, or ',' for free format.
        if (buff[i] == '\n' || buff[i] == ',')
        {
            i++;
            break;
        }
        i++;
    }
    buff[i] = '\0';

    // Continue on new line if end with + and next start with + or *
    if (buff[0] == '+' && buff[i-1] == '\n' &&
        (is.peek() == '+' || is.peek() == '*'))
    {
        if (format == FORMAT::FREE)
        {
            char c;
            while (is.get(c) && c != ',') {}
        }
        else is.readRaw(buff, 8);    // Ignore 1st column

        // Do it again.
        return getColumn(is, size);
    }

    // In free format remove ','
    if (format == FORMAT::FREE)
    {
        // Quirky..
        if (buff[i-1] == ',') buff[i-1] = '\0';
        if (buff[i-1] == '\n') is.putback('\n');
    }

    string data = string(buff);
    // Remove whitespaces
    data.erase(std::remove(data.begin(), data.end(), ' '),
                 data.end());

    return data;
}

string& getEntry(IFstream& is)
{
    // Skip comments
    while (is.good() && is.peek() == '$') is.getLine(nullptr);

    switch (format)
    {
    case FORMAT::SMALL:
    case FORMAT::LARGE:
        entryBuff = getColumn(is, 8);
        break;
    case FORMAT::FREE:
        entryBuff = getColumn(is, 0);
        break;
    }

    if (entryBuff.ends_with('*')) // We have multiline, ignore *.
        entryBuff.removeEnd('*');

    return entryBuff;
}

label getLabel(IFstream& is)
{
    return readLabel(getColumn(is));
}

// Scientific notation sucks in nastran...
scalar getScalar(IFstream& is)
{
    string data = getColumn(is);
    // This is stupid. And probably costly.
    label id = data.find('+', 1); // contains e+
    if (id == -1) id = data.find('-', 1); // contains e-
    if (id != -1) data = data.substr(0, id) + "e" + data.substr(id);

    return readFloat(data);
}

point getPoint(IFstream& is)
{
    // returnt point(get, get, get) will reverse the components. BUG? Or what?
    point pt;
    pt[0] = getScalar(is);
    pt[1] = getScalar(is);
    pt[2] = getScalar(is);
    return pt;
}

// GRID card format, where CP is ignored:
// GRID   ID   CP   X  Y  Z  ...
void readPoints(IFstream& is, DynamicList<point> &points)
{
    label pointi;
    while (true)
    {
        pointi = getLabel(is);
        if (pointi != points.size() + 1)
        {
            FatalErrorInFunction
                << "Points are not in order."
                << exit(FatalError);
        }

        // Ignore CP, we should skip the pointless conversion...
        getLabel(is);
        // Get the 3 coordinate
        points.append(getPoint(is));
        // Ignore the other entries
        skipLine(is); //is.getLine(nullptr);

        if (getEntry(is).starts_with("GRID"))
        {
            continue;
        }
        else
        {
            break;
        }
    }
}

template<cellModel::modelType TYPE>
void readCell
(
    string name,
    IFstream& is,
    DynamicList<cellShape>& cellVerts,
    DynamicList<label>& cellPropIDs
)
{
    const cellModel& model = cellModel::ref(TYPE);
    while (entryBuff == name)
    {
        // Ignore cell ID
        getLabel(is);
        cellPropIDs.append(getLabel(is));
        labelList verts(model.nPoints());
        forAll(verts, i)
        {
            // Vertex indexing starts from 1, but we store from 0.
            verts[i] = getLabel(is) - 1;
        }
        cellVerts.append(cellShape(model, verts, true));

        skipLine(is);
        getEntry(is);
    }
}

template<label nPoints>
void readFaces
(
    string name,
    IFstream& is,
    DynamicList<face>& boundaryFaces,
    DynamicList<label>& facePropIDs
)
{
    while (entryBuff == name)
    {
        getLabel(is); // ignore ID
        facePropIDs.append(getLabel(is));
        face fVerts(nPoints);
        forAll(fVerts, i)
        {
            // Vertex indexing starts from 1, but we store from 0.
            fVerts[i] = getLabel(is) - 1;
        }
        boundaryFaces.append(fVerts);

        skipLine(is); // Get the end of the line
        getEntry(is);
    }
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert nastran dat file to OpenFOAM. Only for small format and Meters."
    );
    argList::noParallel();
    argList::addArgument(".dat file");
    argList::addOption(
        "format",
        "word",
        "Input format. small(default), large, free"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    if (args.found("format"))
    {
        word inpForm = args.getOrDefault<word>("format", "small");
        if (inpForm == "small") 
        {
        // Do nothing. Initialised to small.
        }
        else if (inpForm == "large")
        {
            format = FORMAT::LARGE;
        }
        else if (inpForm == "free")
        {
            format = FORMAT::FREE;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown format: " << inpForm << "."
                << exit(FatalError);
        }
    }

    const auto datName = args.get<fileName>(1);
    IFstream inFile(datName);

    if (!inFile.good())
    {
        FatalErrorInFunction
            << "Cannot open file " << datName
            << exit(FatalError);
    }

    if (!findBulk(inFile))
    {
        FatalErrorInFunction
            << "Cannot find \"BEGIN BULK\" entry."
            << exit(FatalError);
    }

    // Points
    DynamicList<point> points;
    // Cell vertices
    DynamicList<cellShape> cellVerts;
    // Cell property IDs
    DynamicList<label> cellPropIDs;
    // Boundary faces
    DynamicList<face> boundaryFaces;
    // Face property IDs
    DynamicList<label> facePropIDs;
    // Cell zoneIDs
    DynamicList<label> cellZoneIDs;
    // Patch IDs
    DynamicList<label> patchIDs;

    // Read the next entry into the buffer.
    getEntry(inFile);

    Info << "Start reading file." << endl;
    while (inFile.good())
    {
        if (entryBuff == "GRID")
        {
            readPoints(inFile, points);
        }
        else if (entryBuff == "CTETRA")
        {
            readCell<cellModel::modelType::TET>(
                "CTETRA", inFile, cellVerts, cellPropIDs
            );
        }
        else if (entryBuff == "CPYRAM")
        {
            readCell<cellModel::modelType::PYR>(
                "CPYRAM", inFile, cellVerts, cellPropIDs
            );
        }
        else if (entryBuff == "CHEXA")
        {
            readCell<cellModel::modelType::HEX>(
                "CHEXA", inFile, cellVerts, cellPropIDs
            );
        }
        else if (entryBuff == "CTRIA3")
        {
            readFaces<3>("CTRIA3", inFile, boundaryFaces, facePropIDs);
        }
        else if (entryBuff == "PSOLID")
        {
            // add to cell zones
            cellZoneIDs.append(getLabel(inFile));
            skipLine(inFile);
            getEntry(inFile);
        }
        else if (entryBuff == "PSHELL")
        {
            // add to patches
            patchIDs.append(getLabel(inFile));
            skipLine(inFile);
            getEntry(inFile);
        }
        else if (entryBuff.starts_with("ENDDATA")) // EOF did some weird stuff.
        {
            Info << "Finished reading file." << endl;
            break;
        }
        else
        {
            FatalErrorInFunction
                << "Cannot process keyword: \"" << entryBuff.c_str()
                << "\", on line " << inFile.lineNumber() << "."
                << exit(FatalError);
        }
    }

    Info << "Constructing patches. " << endl;
    wordList patchNames;
    forAll(patchIDs, i)
    {
        patchNames.append("patch" + std::to_string(i));
    }

    DynamicList<faceList> patchFaces;
    patchFaces.resize(patchIDs.size());
    forAll(boundaryFaces, i)
    {
        label pI = patchIDs.find(facePropIDs[i]);
        patchFaces[pI].append(boundaryFaces[i]);
    }

    Info << "Constructing the mesh." << endl;
    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        pointField(points),
        cellVerts,
        patchFaces,
        patchNames,
        wordList(patchNames.size(), polyPatch::typeName),
        "defaultFaces",
        polyPatch::typeName,
        wordList()
    );

    
    if (cellZoneIDs.size())
    {
        Info << "Adding cell zones." << endl;

        List<cellZone*> cZones(cellZoneIDs.size());
        
        DynamicList<labelList> cellZones;
        cellZones.resize(cellZoneIDs.size());
        forAll(cellPropIDs, i)
        {
            label zI = cellZoneIDs.find(cellPropIDs[i]);
            cellZones[zI].append(i);
        }
        
        forAll(cellZoneIDs, i)
        {
            cZones[i] = new cellZone
            (
                "cellZone" + std::to_string(i),
                cellZones[i],
                i,
                mesh.cellZones()
            );
        }

        mesh.addZones(List<pointZone*>(), List<faceZone*>(), cZones);
    }

    mesh.removeFiles();
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
