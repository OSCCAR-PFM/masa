/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "solidParticle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidParticle::solidParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            d_ = readScalar(is);   
            is >> U_;
            T_ = readScalar(is); //T
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&d_),
                sizeof(d_) + sizeof(U_) + sizeof(T_) ///T
            );
        }
    }

    // Check state of Istream
    is.check("solidParticle::solidParticle(Istream&)");
}


void Foam::solidParticle::readFields(Cloud<solidParticle>& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);

    IOField<scalar> d(c.fieldIOobject("d", IOobject::MUST_READ));
    c.checkFieldIOobject(c, d);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    IOField<scalar> T(c.fieldIOobject("T", IOobject::MUST_READ));  ///Temp
    c.checkFieldIOobject(c, T);

    label i = 0;
    forAllIter(Cloud<solidParticle>, c, iter)
    {
        solidParticle& p = iter();

        p.d_ = d[i];
        p.U_ = U[i];
        p.T_ = T[i]; //T
        i++;
    }
}


void Foam::solidParticle::writeFields(const Cloud<solidParticle>& c)
{
    particle::writeFields(c);

    label np = c.size();

    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> T(c.fieldIOobject("T", IOobject::NO_READ), np); //T

    label i = 0;
    forAllConstIter(Cloud<solidParticle>, c, iter)
    {
        const solidParticle& p = iter();

        d[i] = p.d_;
        U[i] = p.U_;
        T[i] = p.T_; //T
        i++;
    }

    d.write();
    U.write();
    T.write(); //T
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const solidParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.d_
            << token::SPACE << p.U_
            << token::SPACE << p.T_; ///T
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.d_),
            sizeof(p.d_) + sizeof(p.U_) + sizeof(p.T_) ///
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const solidParticle&)");

    return os;
}


// ************************************************************************* //
