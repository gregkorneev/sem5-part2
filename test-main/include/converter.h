#pragma once
#include <map>
#include <utility>

// Перечисление для единиц измерения температуры
enum TemperatureUnit {
    CELSIUS,
    FAHRENHEIT,
    KELVIN
};

// Перечисление для единиц измерения длины (расширено)
enum LengthUnit {
    METERS,          // метры
    KILOMETERS,      // километры
    MILES,           // мили
    CENTIMETERS,     // сантиметры
    MILLIMETERS,     // миллиметры
    FEET,            // футы
    INCHES,          // дюймы
    YARDS,           // ярды
    NAUTICAL_MILES   // морские мили
};

// Перечисление для валют
enum Currency {
    USD,
    EUR,
    RUB
};

// Объявления функций конвертации
double convertTemperature(double value, TemperatureUnit from, TemperatureUnit to);
double convertLength(double value, LengthUnit from, LengthUnit to);
double convertCurrency(double value, Currency from, Currency to);

// Объявление функции взаимодействия с пользователем
void runConverter();
