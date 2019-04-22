#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 23:16:54 2019

@author: esoh
"""

#firstName = input("First Name: ")
#lastName = input("Last Name: ")
#lastName = lastName.capitalize()
#print("Hello " + firstName + " " + lastName)

#age = 27
#print(age)

#width = 20
#height = 10
#area = width * height /2
#print("The area of the triangle would be %.2f" % area)

#BMI Calculator
#height = float(input("Please enter your height: "))
#weight = float(input("Please provide your weight: "))
#bmi = weight / height**2
#status = input("Are you a heavy weight lifter? [yes/no]: ")
#status = status.lower()
#print("Your BMI is %d kg/m^2" % bmi)
#print("\nHello " + \
#      firstName + \
#      " " + \
#      lastName.upper() + \
#      "\n\nYou are {0:.2f}m tall. You weigh {1:.2f}Kg, and your BMI is {2:.2f}kg/m^2".format(height,weight,bmi))
#
#print(" ")
#if status == "yes":
#    if bmi < 18:
#        print("You are underweight. Please up your diet or see a nutritionist.")
#    elif 18 < bmi <= 25:
#        print("Awesome! You are just fine! Keep up the good diet.")
#    elif 25 < bmi < 30:
#        print("You are overweight. However no fears for a weight lifter like you.")
#    elif bmi >= 30:
#        print("You are obese! Well no worries cos you are a weight lifter :).")
#elif status == "no":
#    if bmi < 18:
#        print("You are underweight. Please up your diet or see a nutritionist.")
#    elif 18 < bmi <= 25:
#        print("Awesome! You are just fine! Keep up the good diet.")
#    elif 25 < bmi < 30:
#        print("You are overweight. Please watch your diet and consider exercises.")
#    elif bmi >= 30:
#        print("You are obese! Please visit a nutritionist and consider exercises!")
## Loan Calculator
#L = float(input("Please enter the cost of your loan: "))
#i = float(input("What is your interest rate? "))
#n = float(input("How many years do you have for your loan? "))
#Monthly = float(L*(i*(1+i)*n)) / float((1+i)*(n-1))
#print("Hello " + firstName + " " + lastName + 
#      " You should pay ${0:.2f} monthly if you are to fully repay your loan of ${1:.2f} on an interest rate of {2:.2f}".format(Monthly,L,i))

# Date and Time
import datetime
currentDate = datetime.date.today()
#print(currentDate)
#print(currentDate.year)
#print(currentDate.month)
#print(currentDate.day)
#print(currentDate.strftime('%d %b,%Y'))
#print(currentDate.strftime('%a %d %B, %Y'))
#print("Could you attend our event " + currentDate.strftime('%A, %B %d') + " in the year " + currentDate.strftime('%Y') + ' "')

#nextBirthday = datetime.datetime.strptime('10/22/2019', '%m/%d/%Y').date()
#print(nextBirthday)
#difference = nextBirthday - currentDate
#print(difference.days)

#deadLine = input("What is the deadline for your project? [mm/dd/yyyy]: ")
#dL = datetime.datetime.strptime(deadLine, '%m/%d/%Y').date()
#difference = dL - currentDate
#weeks = difference.days / 7
#days = difference.days % 7  #This line returns the remainder of the operation (Modulo)
#print("Hey " + firstName + " " + lastName +
#      "! You have " + "{0:.0f}".format(weeks) + "weeks and " + 
#      str(days) + "days to complete your project!")

## The if statement

#freeToaster=False
#deposit = int(input("how much do you want to deposit? "))
#if deposit > 100:
#    freeToaster=True
#
#if freeToaster:
#    print("enjoy your toaster!")
#else:
#    print("Have a nice day!")

#if statements with 'and' and 'or'
#wings = True
#longBeak = True
#if wings and longBeak:
#    print("Darwin's Finch")
#else:
#    print("I don't know this bird")

#het1 = float(input("Please input the heterozygosity for pop1: "))
#het2 = float(input("Please input the heterozygosity for pop2: "))
#het3 = float(input("Please input the heterozygosity for pop3: "))
#het4 = float(input("Please input the heterozygosity for pop4: "))
#het5 = float(input("Please input the heterozygosity for pop5: "))
#het6 = float(input("Please input the heterozygosity for pop6: "))
#withinHet = float(input("Please input the average heterozygosity for pop1 and pop1: "))
#totalHet = float(het1 + het2)
#
#Fst = (totalHet - withinHet) / totalHet
#
#if Fst <= 0.05:
#    print("\nFst = {0:.2f}: Little genetic differenciation".format(Fst))
#elif 0.05 < Fst <= 0.15:
#    print("\nFst = {0:.2f}: Moderate genetic differentiation".format(Fst))
#elif 0.15 < Fst <= 0.25:
#    print("\nFst = {0:.2f}: Great genetic differentiation".format(Fst))
#elif Fst > 0.25:
#    print("\nFst = {0:.2f}: Very great genetic differentiation".format(Fst))

#print(currentDate.strftime('%b'))
#month = str(currentDate.strftime('%b'))
#if month == "Sep" or month == "Apr" \
#    or month == "Jun" or month == "Nov":
#    print("There are 30 days in this month")
#elif month == "Feb":
#    print("There are 28 days in this month or 29 in a leap year")
#else:
#    print("There are 31 days in this month")

# Repeating actions
#import turtle
#turtle.forward(100)
#turtle.right(90)
#turtle.forward(100)
#turtle.right(90)
#turtle.forward(100)
#turtle.right(90)
#turtle.forward(100)
#turtle.right(90)

#times = int(input("Please input the number of steps you want turtle to run: "))
#
#for steps in range(times):
#    turtle.forward(100)
#    turtle.right(360/times)
#    for moresteps in range(times):
#        turtle.forward(50)
#        turtle.right(360/times)

#for steps in range(1,4):
#    print(steps)
    
#for steps in [1,2,3,4,5]:
#    print(steps)
#
#for color in ['red', 'green', 'blue', 'brown']:
#    turtle.color(color)
#    turtle.forward(100)
#    turtle.left(90)

#answer = "0"
#while answer != "4":
#    answer = input("What is 2 + 2? ")
#    print("\nSorry that is the wrong answer. Please try again")
#print("Yes 2 + 2 is 4")

#counter = 0
#while counter < 4:
#    turtle.forward(100)
#    turtle.right(90)
#    counter = counter+1
#    counter += 1

#Reading values from files (lists)
#guest = ["Kev", "Jay", "Kay", "Jo"]
#guest[0] = 'Rev'
#guest.append('Steve')
#guest.remove('Jo')
#del guest[-2]
#guest.append('Colin')
#print(guest.index('Jay'))
#print(guest[int(input("Which index do u want to print? "))])

#for currentGuest in guest:
#    print(currentGuest)

#for steps in range(len(guest)):
#    print(guest[steps])
#    
#guest.sort()

#guests = [ ]
#name = " "
#while name != "Done":
#    name = input("Enter name of invitee (enter DONE if no more names): ").capitalize()
#    guests.append(name)
#
#guests.remove('Done')
#guests.sort()
#for guest in guests:
#    print(guest)
    
# Writing out a list to a file
#myFile1 = "GuestList.txt"
#READ = "r" 
#WRITE = "w"
#APPEND = "a"
#READWRITE = "w+"
#
#myFile1 = open(myFile1, WRITE)
#myFile1.write("Hey there\n")
#myFile1.write("How are you doing?\n")
#myFile1.write("I love you\n")
#myFile1.close()
#
#myFile2 = "GuestList.csv"
#myFile2 = open(myFile2, READWRITE)
#myFile2.write("Kev, 1\n")
#myFile2.write("Kum, 3\n")
#myFile2.write("Jay, 4")
#myFile2.close()
#
#data = input("Please input file info: ")
#myFile3 = open('inputData.txt', READWRITE)
#myFile3.write(data)
#myFile3.close()

#Read from a file
# First open the file
#animalFile = open("inputData.txt", "r")
# Then read the contents of the file
#firstAnimal = animalFile.readline()
#print(firstAnimal)
#secondAnimal = animalFile.readline()
#print(secondAnimal)

#import csv
#fileName = "inputData.txt"
#accessMode = "r"
#with open(fileName, accessMode) as myCSVFile:
#    dataFromFile = csv.reader(myCSVFile)
#    
#    for currentRow in dataFromFile :
#        print(', ' .join(currentRow))
##        for currentWord in currentRow:
##            print(currentWord)
#
#newFile = "newChallenge.txt"
#READWRITE = "r"
#newChallenge = open(newFile, accessMode)
#newChallenge.write("Susan,32\nMike,33\nChristopher,80\nBill Gates,100")
#newChallenge.close()

#with open("newChallenge.txt", accessMode) as challengeFile:
#    challengeData = csv.reader(challengeFile)
#    for currentRow in challengeData:
#        print(', ' .join(currentRow))
        
#Creating Functions
#Define a first function that will contain all the functions you will create
#def main():
#    printMessage()
#    printMyName()
#    greeting = "Hello!"
#    printGreeting(greeting, input("Please enter your name: ").capitalize())
#    return
#
#def printMessage():
#    print("Hello World!")
#    return
#
#def printMyName():
#    listOfCities = ["Toronto\n", "New York\n", "Quebec\n", "Montreal\n", "Paris\n"]
#    myName = input("Please enter your name: ").capitalize()
#    print("\nHello " + myName + "! Welcome")
#    print("Which of these cities do you prefer?" + "\n")
#    print('' .join(listOfCities))
#    return
#
#def printGreeting(greeting, name):
#    message = name + ', ' + greeting
#    print(message)
#    return 
#
#printGreeting(input("Please enter your message: "), input("Please provide your name: "))
#main()

#import pandas as pd
#import numpy as np

from Bio import Entrez
import datetime
import csv
Entrez.email = "kevin.esoh@students.jkuat.ac.ke"
apiKey=''.join(open("api_key.txt", "r"))
api_key = apiKey.lower()
currentDate = datetime.date.today()
currentTime = datetime.datetime.now()

def main():
    dataBase = getDataBase()
    searchTerm = getSearchTerm()
    dbSearch(dataBase, searchTerm)
    with open("idList.txt", "r+") as newFile:
        newIdList = csv.reader(newFile)
        for eachID in newIdList:
            idList = ','.join(eachID)
    dbSummary(dataBase, idList)
    return

def getDataBase():
    dataBase = input("Please provide the NCBI database you wish to search: ").lower()
    return dataBase

def getSearchTerm():
    searchTerm = input("Enter you search word(s): ").lower()
    return searchTerm

def dbSearch(dataBase, searchTerm):
    handle = Entrez.esearch(db=dataBase, retmax=30, term=searchTerm, retmode="uid")
    record = Entrez.read(handle)
    #record2 = Entrez.read(handle, validate=True)
    handle.close()
    print(" ")
    print("Your query: " + record["QueryTranslation"])
    print(" ")
    print("Total count: " + record["Count"] + \
          " Accessed: " + currentDate.strftime('%a %d %B, %Y') + \
          " At: " + currentTime.strftime('%H:%M:%S'))
    idList = record["IdList"]
    print(" ")
    print(idList)
    newList = open("idList.txt", "w+")
    newList.write(',' .join(idList))
    newList.close()
    #print(idList)
    #print(',' .join(idList))
    #print(" ")
    #print(record["Title"])
    return

def dbSummary(dataBase, idList=[]):
    handle = Entrez.esummary(db=dataBase, id=idList, rettype="xml")
    #records = Entrez.read(handle)
    records = Entrez.parse(handle)
    for record in records:
        print(" ")
        #print(record)
        print(record['Id'] + " | " + \
                record['Title'] + " | " + \
                record['AuthorList'][0] + " | " + \
                record['PubDate'])
    handle.close()
    print(" ")
    return
main()
