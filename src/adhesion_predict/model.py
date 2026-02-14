"""Model training and prediction functions."""

import pickle
import sys

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, train_test_split


def train_classifier(x, y, test_size=0.2, random_state=42):
    """Train logistic regression classifier on embeddings.

    Args:
        x: Feature matrix.
        y: Labels.
        test_size: Proportion for test set.
        random_state: Random seed.

    Returns:
        Tuple of (trained classifier, test accuracy).
    """
    x_train, x_test, y_train, y_test = train_test_split(
        x, y, test_size=test_size, random_state=random_state, stratify=y
    )

    print(f"Training set size: {len(x_train)}")
    print(f"Test set size: {len(x_test)}")

    print("Training classifier...")
    clf = LogisticRegression(max_iter=1000, random_state=random_state)
    clf.fit(x_train, y_train)

    train_accuracy = clf.score(x_train, y_train)
    test_accuracy = clf.score(x_test, y_test)

    print(f"  Training accuracy: {train_accuracy:.3f}")
    print(f"  Test accuracy: {test_accuracy:.3f}")

    cv_scores = cross_val_score(clf, x_train, y_train, cv=5, scoring="accuracy")
    print(f"  Cross-validation accuracy: {cv_scores.mean():.3f} (+/- {cv_scores.std() * 2:.3f})")

    print("Training final model on all data...")
    clf_final = LogisticRegression(max_iter=1000, random_state=random_state)
    clf_final.fit(x, y)

    return clf_final, test_accuracy


def save_model(classifier, output_path):
    """Save trained classifier to file.

    Args:
        classifier: Trained classifier.
        output_path: Path to save the model.
    """
    with open(output_path, "wb") as f:
        pickle.dump(classifier, f)
    print(f"Model saved to {output_path}")


def load_model(model_path):
    """Load trained classifier from file.

    Args:
        model_path: Path to the saved model.

    Returns:
        Loaded classifier.
    """
    try:
        with open(model_path, "rb") as f:
            classifier = pickle.load(f)
        return classifier
    except FileNotFoundError:
        print(f"Error: Model file not found at {model_path}", file=sys.stderr)
        sys.exit(1)


def predict(classifier, X):
    """Make predictions using trained classifier.

    Args:
        classifier: Trained classifier.
        X: Feature matrix.

    Returns:
        Predictions array.
    """
    return classifier.predict(X)


def predict_proba(classifier, X):
    """Get prediction probabilities.

    Args:
        classifier: Trained classifier.
        X: Feature matrix.

    Returns:
        Array of prediction probabilities for each class.
    """
    return classifier.predict_proba(X)
